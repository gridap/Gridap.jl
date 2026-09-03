# RotatingPLambda/BarycentricPΛTrimmedBases.jl
#
# Trimmed P_r⁻Λ¹ basis: same architecture as
# RotatingPLambda/BarycentricPΛBases.jl, but built from the Whitney
# (directional) 1-forms
#
#   ϕ(ξ;e1,e2) = ξ_e1 dξ_e2 − ξ_e2 dξ_e1
#
# instead of the rotating ψ. Spanning set/basis: f a face of dim(f)≥1,
# e=(e1,e2) a pair of distinct vertices of f, |α|=r-1, [s(α)]∪e=f, basis
# filter e1<e2 and e1=min(f).
#
# Unlike ψ, ϕ is NOT a constant-coefficient directional form: ϕ(ξ;e1,e2) is
# itself linear in ξ. So the basis function
#   w(ξ;F,e,α) = BB(ξ;α)·ϕ(ξ;e1,e2) = BB(ξ;α)·ξ_e1·dξ_e2 − BB(ξ;α)·ξ_e2·dξ_e1
# is evaluated as (scalar Bernstein value)·(barycentric coordinate)·(constant
# physical dx-form), combining two such terms, rather than a single constant
# direction times a scalar Bernstein value as in the untrimmed case.
#
# Only K=1 forms are supported.

##############################
# Bubble indices (per face F) #
##############################

# One basis function of face F: (w, e=(e1,e2), α, α_id)
const TrimmedBubbleFunction = Tuple{Int,Tuple{Int,Int},Vector{Int},Int}

# One face's bubble space: (F, [bubble functions])
const TrimmedBubble = Tuple{Vector{Int},Vector{TrimmedBubbleFunction}}

"""
    _trimmed_F_bubble_functions(r, D, F, w)

Bubble functions of face `F` (a sorted vertex set of the D-simplex): all
`((e1,e2),α)` with `e1 = min(F)` fixed (the basis filter), `e2 ∈ F∖{e1}`,
`|α|=r-1`, and `supp(α)∪{e1,e2}=F`. `w` is the running basis-function index.
"""
function _trimmed_F_bubble_functions(r::Int, D::Int, F::Vector{Int}, w::Int)
  bubble_functions = TrimmedBubbleFunction[]
  e1 = minimum(F)
  for e2 in F
    e2 == e1 && continue
    for α in bernstein_terms(r-1, D)
      supp = [i for i in 1:(D+1) if α[i] > 0]
      issetequal(union(supp, (e1,e2)), F) || continue
      w += 1
      push!(bubble_functions, (w, (e1,e2), α, bernstein_term_id(α)))
    end
  end
  bubble_functions
end

"""
    trimmed_PΛ_bubbles(r, D)

Generates the bubble indices of the trimmed P_r⁻Λ¹ basis, grouped by face,
exactly like [`rotating_PΛ_bubbles`](@ref) does for the untrimmed case.
"""
function trimmed_PΛ_bubbles(r::Int, D::Int)
  bubbles = TrimmedBubble[]
  w = 0
  for d in 1:D   # dim(f)≥1 always: a pair (e1,e2) needs ≥2 vertices in f
    for F in combinations(1:(D+1), d+1)
      bf = _trimmed_F_bubble_functions(r, D, F, w)
      isempty(bf) && continue
      push!(bubbles, (F, bf))
      w += length(bf)
    end
  end
  bubbles
end

###################################
# ϕ(e1,e2) ambient coefficient    #
###################################

"""
    _trimmed_ambient_phi(e1, e2, N, λ)

ϕ(ξ;e1,e2) = λ_e1 dλ_e2 − λ_e2 dλ_e1, evaluated at the ambient barycentric
point/symbol vector `λ` (length `N`), as a coefficient vector in the ambient
frame dλ₁,…,dλ_N. Unlike [`_rotating_ambient_psi`](@ref), this genuinely
depends on the evaluation point — `λ` here may be numeric (for evaluation) or
symbolic (for [`print_forms`](@ref)).
"""
function _trimmed_ambient_phi(e1::Int, e2::Int, N::Int, λ)
  c = zeros(eltype(λ), N)
  c[e2] += λ[e1]
  c[e1] -= λ[e2]
  c
end

##########################################
# TrimmedPΛBasis{D,V,C,K} nD polynomial  #
##########################################

"""
    TrimmedPΛBasis{D,V,C,K} <: PolynomialBasis{D,V,Bernstein}

Trimmed P_r⁻Λ¹ basis on the D-simplex: `B_{α,r-1}(λ) · ϕ(e1,e2)`, with `ϕ` the
Whitney 1-form (see comments above).

- `V = DifferentialFormValue{1,D,T,D}`,
- `C` the number of basis polynomials,
- `K` the polynomial order of the underlying scalar Bernstein basis (= r−1).
"""
struct TrimmedPΛBasis{D,V,C,K} <: PolynomialBasis{D,V,Bernstein}
  scalar_bernstein_basis :: BernsteinBasisOnSimplex{D,Float64,K}  # degree r-1
  e1s     :: Vector{Int}             # e1, indexed by w
  e2s     :: Vector{Int}             # e2, indexed by w
  Je1     :: Vector{V}               # physical dx-form for unit dλ_e1, indexed by w
  Je2     :: Vector{V}               # physical dx-form for unit dλ_e2, indexed by w
  α_ids   :: Vector{Int}             # bernstein_term_id(α), indexed by w
  ∇λ      :: Vector{VectorValue{D,Float64}}  # ∂λ_i/∂x (constant), indexed by vertex i
  nrm     :: Vector{Float64}         # 1/multinomial(α): BARE monomial λ^α, not B_α.
                                     # The ±1 two-term rotation law (shift ρ changes
                                     # the multiset of α) holds for bare monomials
                                     # only; Bernstein normalisation would scale the
                                     # hit entries by multinomial ratios.
  bubbles :: Vector{TrimmedBubble}   # grouped by face, for inspection / printing

  function TrimmedPΛBasis{D}(::Type{T}, r::Int, vertices=nothing) where {D,T}
    @assert T<:Real "T needs to be <:Real, got $T"
    @assert r ≥ 1 "r must be ≥ 1, got $r (|α|=r-1 requires r≥1)"

    N = D + 1
    bubbles = trimmed_PΛ_bubbles(r, D)
    C = isempty(bubbles) ? 0 : bubbles[end][2][end][1]

    L = D # binomial(D,1)
    V = DifferentialFormValue{1,D,T,L}

    b = BernsteinBasisOnSimplex{D}(Float64, r-1, vertices)
    K = get_order(b)

    e1s   = Vector{Int}(undef, C)
    e2s   = Vector{Int}(undef, C)
    Je1   = Vector{V}(undef, C)
    Je2   = Vector{V}(undef, C)
    α_ids = Vector{Int}(undef, C)
    nrm   = Vector{Float64}(undef, C)
    J     = _rotating_lambda_jacobian(Val(D), b.x_to_λ)
    ∇λ    = [VectorValue{D,Float64}(NTuple{D,Float64}(J[:,i])) for i in 1:N]

    for (F, bubble_functions) in bubbles
      for (w, e, α, α_id) in bubble_functions
        e1, e2   = e
        e1s[w]   = e1
        e2s[w]   = e2
        Je1[w]   = DifferentialFormValue{1,D}(NTuple{L,T}(T.(J[:,e1])))
        Je2[w]   = DifferentialFormValue{1,D}(NTuple{L,T}(T.(J[:,e2])))
        α_ids[w] = α_id
        nrm[w]   = 1.0 / multinomial(α...)
      end
    end

    new{D,V,C,K}(b, e1s, e2s, Je1, Je2, α_ids, ∇λ, nrm, bubbles)
  end
end

"""
    TrimmedPΛBasis(::Val{D}, T, r, vertices=nothing)

Constructor for [`TrimmedPΛBasis`](@ref).
"""
TrimmedPΛBasis(::Val{D},::Type{T},r,vertices=nothing) where {D,T} =
  TrimmedPΛBasis{D}(T,r,vertices)

get_bubbles(b::TrimmedPΛBasis) = b.bubbles
get_order(b::TrimmedPΛBasis) = get_order(b.scalar_bernstein_basis) + 1
get_orders(b::TrimmedPΛBasis{D}) where D = ntuple(_ -> get_order(b), D)

Base.size(::TrimmedPΛBasis{D,V,C}) where {D,V,C} = (C, )

# ── evaluate! / return_cache: PolynomialBasis low-level API ────────────────

# Cache layout mirrors _BaryPΛBasis (see RotatingPLambda/BarycentricPΛBases.jl).
function _return_cache(b::TrimmedPΛBasis, x, ::Type{G}, ::Val{N_deriv}) where {G,N_deriv}
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

function _setsize!(b::TrimmedPΛBasis, np, ω, t...)
  ndof = length(b)
  ndof_bernstein = length(b.scalar_bernstein_basis)
  setsize!(ω,(np,ndof))
  setsize!(t[1],(ndof_bernstein,))
  if length(t) > 1
    setsize!(t[end],(1,ndof_bernstein))
  end
end

_get_parameters(b::TrimmedPΛBasis) = Val(get_order(b.scalar_bernstein_basis))

function _evaluate_nd!(
  b::TrimmedPΛBasis{D}, x,
  ω::AbstractMatrix, i, cB,
  ::Val{r}) where {D,r}

  λ = _cart_to_bary(x, b.scalar_bernstein_basis.x_to_λ)

  cB[1] = 1
  _downwards_de_Casteljau_nD!(cB,λ,Val(r),Val(D))

  @inbounds for w in 1:length(b)
    bb = b.nrm[w] * cB[b.α_ids[w]]   # bare monomial λ^α = B_α / multinomial(α)
    ω[i,w] = bb*λ[b.e1s[w]]*b.Je2[w] - bb*λ[b.e2s[w]]*b.Je1[w]
  end
end

# Product rule for w = nrm·B_α·(λ_{e1}·Je2 − λ_{e2}·Je1) with constant Je:
#   ∇w = nrm·[ (λ_{e1}∇B_α + B_α∇λ_{e1}) ⊗ Je2 − (λ_{e2}∇B_α + B_α∇λ_{e2}) ⊗ Je1 ]
function _gradient_nd!(
  b::TrimmedPΛBasis{D}, x,
  ∇ω::AbstractMatrix{G}, i, cB,
  ∇B::AbstractMatrix{<:VectorValue{D}},
  s::MVector{D},
  ::Val{r}) where {D,G,r}

  # gradients of all scalar Bernstein polynomials (fills ∇B row 1; cB is used
  # as scratch at degree r−1 internally, so recompute the values afterwards)
  _gradient_nd!(b.scalar_bernstein_basis, x, ∇B, 1, cB, nothing, s, Val(r))

  λ = _cart_to_bary(x, b.scalar_bernstein_basis.x_to_λ)
  cB[1] = 1
  _downwards_de_Casteljau_nD!(cB,λ,Val(r),Val(D))

  @inbounds for w in 1:length(b)
    e1, e2 = b.e1s[w], b.e2s[w]
    bb  = b.nrm[w] * cB[b.α_ids[w]]
    ∇bb = b.nrm[w] * ∇B[1,b.α_ids[w]]
    g1  = λ[e1]*∇bb + bb*b.∇λ[e1]
    g2  = λ[e2]*∇bb + bb*b.∇λ[e2]
    ∇ω[i,w] = outer(g1, b.Je2[w]) - outer(g2, b.Je1[w])
  end
end

# ── Pretty table ─────────────────────────────────────────────────────────

"""
    print_indices(b::TrimmedPΛBasis, out::IO=stdout)

Pretty table of the trimmed P_r⁻Λ¹ basis, one row per basis function `w`:
face `F`, vertex pair `e=(e1,e2)`, multi-index `α`, the bare Bernstein
monomial `λ^α`, and the physical Jacobian columns `Je1,Je2` (the Cartesian
`dx`-form for unit `dλ_e1`/`dλ_e2`). Mirrors [`print_indices`](@ref).
"""
function print_indices(b::TrimmedPΛBasis{D}, out::IO=stdout) where D
  println(out, "TrimmedPΛBasis{D=$D, r=$(get_order(b))}: dim = $(length(b))")
  println(out,
    rpad("w",4), rpad("F",10), rpad("e",7), rpad("α",12),
    rpad("B_α(λ)",18), "Je1, Je2 (dx¹,…,dx^$D)")
  for (F, bubble_functions) in b.bubbles
    for (w, e, α, _) in bubble_functions
      mono = _monomial_string(α)   # bare monomial
      println(out,
        rpad("$w",4), rpad(join(F,","),10), rpad("$e",7), rpad("$(Tuple(α))",12),
        rpad(mono,18), "$(b.Je1[w]), $(b.Je2[w])")
    end
  end
end

Base.show(io::IO, b::TrimmedPΛBasis) = print_indices(b, io)
