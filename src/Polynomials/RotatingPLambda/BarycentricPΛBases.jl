# RotatingPLambda/BarycentricPΛBases.jl
#
# Rotating P_rΛ¹ basis: same architecture as `BarycentricPΛBasis`
# (Bernstein-Bézier scalar basis + precomputed per-basis-function coefficient
# table), but built from the directional 1-form
#
#   ψ(f,k,α) = dλᵏ − (𝟙[k∈supp(α)]/|supp(α)|) · Σ_{i∈f} dλⁱ
#
# instead of the AFW geometric decomposition. Its basis functions evaluate to
# genuine `DifferentialFormValue{1,D}` 1-forms, and its behaviour under vertex
# relabellings admits a closed-form change of basis (see
# RotatingPLambda/PΛRotations.jl), which is what makes conformity on
# arbitrarily ordered simplicial meshes cheap.
#
# Only K=1 forms are supported: the formula above picks out a single vertex k.

##############################
# Bubble indices (per face F) #
##############################

# One basis function of face F: (w, k, α, α_id)
const RotatingBubbleFunction = Tuple{Int,Int,Vector{Int},Int}

# One face's bubble space: (F, [bubble functions])
const RotatingBubble = Tuple{Vector{Int},Vector{RotatingBubbleFunction}}

"""
    _rotating_F_bubble_functions(r, D, F, w)

Bubble functions of face `F` (a sorted vertex set of the D-simplex): all
`(k,α)` with `|α|=r`, `k∈F`, `supp(α)∪{k}=F`, and the basis filter
`αᵢ=0` for `i<min(F\\k)`. `w` is the running basis-function index.
"""
function _rotating_F_bubble_functions(r::Int, D::Int, F::Vector{Int}, w::Int)
  bubble_functions = RotatingBubbleFunction[]
  for α in bernstein_terms(r, D)
    supp = [i for i in 1:(D+1) if α[i] > 0]
    for k in F
      issetequal(union(supp, (k,)), F) || continue
      rest = setdiff(F, (k,))
      j = isempty(rest) ? 0 : minimum(rest) - 1
      all(==(0), α[1:j]) || continue
      w += 1
      push!(bubble_functions, (w, k, α, bernstein_term_id(α)))
    end
  end
  bubble_functions
end

"""
    rotating_PΛ_bubbles(r, D)

Generates the bubble indices of the rotating P_rΛ¹ basis, grouped by face,
exactly like [`PΛ_bubbles`](@ref) does for the AFW basis:

    for (F, bubble_functions) in rotating_PΛ_bubbles(r,D)
      for (w, k, α, α_id) in bubble_functions
        # ...
      end
    end
"""
function rotating_PΛ_bubbles(r::Int, D::Int)
  bubbles = RotatingBubble[]
  w = 0
  for d in 1:D   # K=1: face dimension ranges over 1:D (edges up to the full simplex)
    for F in combinations(1:(D+1), d+1)
      bf = _rotating_F_bubble_functions(r, D, F, w)
      isempty(bf) && continue
      push!(bubbles, (F, bf))
      w += length(bf)
    end
  end
  bubbles
end

###############################
# ψ(f,k,α) ambient coefficient #
###############################

"""
    _rotating_ambient_psi(F, k, α, N)

ψ(F,k,α) = dλᵏ − (𝟙[k∈supp(α)]/|supp(α)|) · Σ_{i∈F} dλⁱ, as a coefficient
vector in the ambient N-dim barycentric frame (dλ¹,…,dλᴺ).
"""
function _rotating_ambient_psi(F::Vector{Int}, k::Int, α::Vector{Int}, N::Int)
  c = zeros(Float64, N)
  c[k] += 1.0
  supp = [i for i in 1:N if α[i] > 0]
  if !isempty(supp)
    w = (k in supp) ? 1.0/length(supp) : 0.0
    for i in F
      c[i] -= w
    end
  end
  c
end

# Jacobian J[a,i] = ∂λⁱ/∂xₐ (D×N), used to contract the ambient ψ coefficients
# down to genuine physical Cartesian dx-form coefficients. It comes from
# BernsteinBasisOnSimplex's own `x_to_λ` (λ = M*[1;x]), exactly like
# _update_φ_αF! does for BarycentricPΛBasis.
function _rotating_lambda_jacobian(::Val{D}, M) where D
  N = D + 1
  J = zeros(Float64, D, N)
  for a in 1:D, i in 1:N
    J[a,i] = M[i,a+1]
  end
  J
end

##########################################
# RotatingPΛBasis{D,V,C,K} nD polynomial  #
##########################################

"""
    RotatingPΛBasis{D,V,C,K} <: PolynomialBasis{D,V,Bernstein}

Rotating P_rΛ¹ basis on the D-simplex: `B_{α,r}(λ) · ψ(F,k,α)`, with `ψ` the
directional 1-form (see comments above).

- `V = DifferentialFormValue{1,D,T,D}`,
- `C` the number of basis polynomials,
- `K` the polynomial order of the underlying scalar Bernstein basis (= r).
"""
struct RotatingPΛBasis{D,V,C,K} <: PolynomialBasis{D,V,Bernstein}
  scalar_bernstein_basis :: BernsteinBasisOnSimplex{D,Float64,K}
  Ψ       :: Vector{V}              # physical 1-form coefficient, indexed by w
  α_ids   :: Vector{Int}            # bernstein_term_id(α), indexed by w
  bubbles :: Vector{RotatingBubble} # grouped by face, for inspection / printing

  function RotatingPΛBasis{D}(::Type{T}, r::Int, vertices=nothing) where {D,T}
    @assert T<:Real "T needs to be <:Real, got $T"
    @assert r ≥ 1 "r must be ≥ 1, got $r (r=0 admits no geometric decomposition for 1-forms)"

    N = D + 1
    bubbles = rotating_PΛ_bubbles(r, D)
    C = isempty(bubbles) ? 0 : bubbles[end][2][end][1]

    L = D # binomial(D,1)
    V = DifferentialFormValue{1,D,T,L}

    b = BernsteinBasisOnSimplex{D}(Float64, r, vertices)
    K = get_order(b)

    Ψ     = Vector{V}(undef, C)
    α_ids = Vector{Int}(undef, C)
    J     = _rotating_lambda_jacobian(Val(D), b.x_to_λ)

    for (F, bubble_functions) in bubbles
      for (w, k, α, α_id) in bubble_functions
        c    = _rotating_ambient_psi(F, k, α, N)
        phys = J * c
        Ψ[w]     = DifferentialFormValue{1,D}(NTuple{L,T}(T.(phys)))
        α_ids[w] = α_id
      end
    end

    new{D,V,C,K}(b, Ψ, α_ids, bubbles)
  end
end

"""
    RotatingPΛBasis(::Val{D}, T, r, vertices=nothing)

Constructor for [`RotatingPΛBasis`](@ref).
"""
RotatingPΛBasis(::Val{D},::Type{T},r,vertices=nothing) where {D,T} =
  RotatingPΛBasis{D}(T,r,vertices)

get_bubbles(b::RotatingPΛBasis) = b.bubbles
get_order(b::RotatingPΛBasis) = get_order(b.scalar_bernstein_basis)
get_orders(b::RotatingPΛBasis{D}) where D = ntuple(_ -> get_order(b), D)

Base.size(::RotatingPΛBasis{D,V,C}) where {D,V,C} = (C, )

# ── evaluate! / return_cache: PolynomialBasis low-level API ────────────────

# Cache layout mirrors _BaryPΛBasis: (r, s, cB) for values,
# (r, s, cB, ∇B) for first derivatives (∇B = gradients of all scalar Bernstein
# polynomials at one point, s = MVector scratch of the scalar kernel).
function _return_cache(b::RotatingPΛBasis, x, ::Type{G}, ::Val{N_deriv}) where {G,N_deriv}
  T  = eltype(G)
  np = length(x)
  ndof = length(b)
  ndof_bernstein = length(b.scalar_bernstein_basis)

  r  = CachedArray(zeros(G,(np,ndof)))
  cB = CachedVector(zeros(T,ndof_bernstein))
  if N_deriv > 0
    DB = T
    xi = testitem(x)
    for _ in 1:N_deriv
      DB = gradient_type(DB,xi)
    end
    t = (( nothing for _ in 2:N_deriv)..., CachedArray(zeros(DB,(1,ndof_bernstein))))
    s = MArray{Tuple{size(DB)...},T}(undef)
  else
    t = ()
    s = nothing
  end
  (r, s, cB, t...)
end

function _setsize!(b::RotatingPΛBasis, np, ω, t...)
  ndof = length(b)
  ndof_bernstein = length(b.scalar_bernstein_basis)
  setsize!(ω,(np,ndof))
  setsize!(t[1],(ndof_bernstein,))
  if length(t) > 1
    setsize!(t[end],(1,ndof_bernstein))
  end
end

_get_parameters(b::RotatingPΛBasis) = Val(get_order(b))

function _evaluate_nd!(
  b::RotatingPΛBasis{D}, x,
  ω::AbstractMatrix, i, cB,
  ::Val{r}) where {D,r}

  λ = _cart_to_bary(x, b.scalar_bernstein_basis.x_to_λ)

  cB[1] = 1
  _downwards_de_Casteljau_nD!(cB,λ,Val(r),Val(D))

  @inbounds for w in 1:length(b)
    ω[i,w] = cB[b.α_ids[w]] * b.Ψ[w]
  end
end

# ∇(B_α·Ψ) = ∇B_α ⊗ Ψ (Ψ is a constant-coefficient 1-form).
function _gradient_nd!(
  b::RotatingPΛBasis{D}, x,
  ∇ω::AbstractMatrix{G}, i, cB,
  ∇B::AbstractMatrix{<:VectorValue{D}},
  s::MVector{D},
  ::Val{r}) where {D,G,r}

  _gradient_nd!(b.scalar_bernstein_basis, x, ∇B, 1, cB, nothing, s, Val(r))

  @inbounds for w in 1:length(b)
    ∇ω[i,w] = outer(∇B[1,b.α_ids[w]], b.Ψ[w])
  end
end

# ── Pretty table ─────────────────────────────────────────────────────────

"""
    print_indices(b::RotatingPΛBasis, out::IO=stdout)

Pretty table of the rotating P_rΛ¹ basis, one row per basis function `w`:
face `F`, vertex `k`, multi-index `α`, the Bernstein monomial `B_α(λ)`, and
the physical `ψ` (the Cartesian `dx`-form coefficients, i.e. `b.Ψ[w]`).
Mirrors [`print_indices`](@ref) for the AFW `BarycentricPΛBasis`.
"""
function print_indices(b::RotatingPΛBasis{D}, out::IO=stdout) where D
  N = D + 1
  println(out, "RotatingPΛBasis{D=$D, r=$(get_order(b))}: dim = $(length(b))")
  println(out,
    rpad("w",4), rpad("F",10), rpad("k",4), rpad("α",12),
    rpad("B_α(λ)",18), "ψ (dx¹,…,dx^$D)")
  for (F, bubble_functions) in b.bubbles
    for (w, k, α, _) in bubble_functions
      mono = _monomial_string(α; coeff=multinomial(α...))
      println(out,
        rpad("$w",4), rpad(join(F,","),10), rpad("$k",4), rpad("$(Tuple(α))",12),
        rpad(mono,18), b.Ψ[w])
    end
  end
end

Base.show(io::IO, b::RotatingPΛBasis) = print_indices(b, io)

"""
    print_forms(b::BarycentricPΛBasis,  out::IO=stdout)
    print_forms(b::BarycentricPmΛBasis, out::IO=stdout)
    print_forms(b::RotatingPΛBasis,     out::IO=stdout)
    print_forms(b::TrimmedPΛBasis,      out::IO=stdout)

Print each basis function of `b` as an ambient barycentric differential form,
in terms of dλ¹,…,dλ^{D+1} — the frame the φ, ψ and ϕ formulas are stated in —
with symbolic polynomial coefficients.

For a `BarycentricPΛBasis` the form is `Bα(λ)` times the wedge of the direction
1-forms [`_update_φ_αF!`](@ref) indexed by `J`, and the polynomial degree must be
at least 1. For a `BarycentricPmΛBasis` it is the Whitney form `φ^J` scaled by
`Bα(λ)` for `flavor=:AFW`, by the bare monomial `λ^α` for `flavor=:BMM`.

Requires the Symbolics package to be loaded: the methods are provided by the
GridapSymbolicsExt package extension.
"""
function print_forms end
