#####################################
# P⁻LambdaBasis nD polynomial bases #
#####################################

"""
struct PmLambdaBasis{D,V,L,B} <: PolynomialBasis{D,V,Bernstein}

Finite Element Exterior Calculus polynomial basis for the spaces P⁻`ᵣ`Λ`ᴷ` on
`D`-dimensional simplices, but with polynomial forms explicitely transformed
into vectors using the standard equivalence with usual vector calculus (involving
the hodge star operator ⋆ and the flat map ♭).

`V` is `VectorValue{L,T}` where
`L` is binomial(`D`,`k`)
`B` is the concrete type of the `BernsteinBasisOnSimplex` necessary for the evaluation of the polynomials.

The number of basis polynomials is binomial(`r`+`k`-1,`k`)*binomial(`D`+`r`,`D`-`k`).

Reference: D. N. Arnold and A. Logg, Periodic Table of the Finite Elements, SIAM News, vol. 47 no. 9, November 2014
"""
struct PmLambdaBasis{D,V,LN,B} <: PolynomialBasis{D,V,Bernstein}
  r::Int
  k::Int
  scalar_bernstein_basis::B
  m::SVector{LN,V}
  _indices::PLambdaIndices

  function PmLambdaBasis{D}(::Type{T}, r, k, vertices=nothing;
      diff_geo_calculus_style=false, indices=nothing) where {D,T}

    @check T<:Real "T needs to be <:Real since represents the scalar type, got $T"
    @check k in 0:D "The form order k must be in 0:D, got $k"
    @check r > 0    "The polynomial order r must be positive, got $r"
    if !isnothing(vertices)
      @check length(vertices) == D+1 "$D+1 vertices are required to define a $D-dim simplex, got $(length(vertices))"
      @check eltype(vertices) <: Point{D} "Vertices should be of type <:Point{$D}, got $(eltype(vertices))"
    end

    indices = _generate_PmΛ_indices(r,k,D,diff_geo_calculus_style,indices)

    L = binomial(D,k) # Number of components of a basis form

    if diff_geo_calculus_style
      @notimplemented
    else
      (D>3 && ( 1 < k < D-1)) && @unreachable "Vector calculus proxy of differential form bases not available for D=$D and k=$k"
      V = VectorValue{L,T}
    end

    b = BernsteinBasisOnSimplex{D}(T, r, vertices)
    B = typeof(b)
    LN = binomial(D+1,k) # Number of k faces of a D dimensional tetrahedron
    m = zero(MVector{LN,V})
    _compute_PmΛ_basis_form_coefficient!(m,k,D,b,vertices,indices)

    if isone(L) && !diff_geo_calculus_style
      V = T
      m = reinterpret(T, m)
    end

    new{D,V,LN,B}(r,k,b,m,indices)
  end
end

function PmLambdaBasis(::Val{D},::Type{T},r,k,vertices=nothing; diff_geo_calculus_style=false) where {D,T}
  PmLambdaBasis{D}(T,r,k,vertices; diff_geo_calculus_style)
end

"""
    _DG_calculus_style(::V) = false

Temporary API to signal that the coefficient of `V`-valued polynomial forms
should be transformed into classic vector calulus components using Hodge star
operator.
"""
_DG_calculus_style(::V) where V = false

get_FEEC_poly_degree(b::PmLambdaBasis) = b.r
get_FEEC_form_degree(b::PmLambdaBasis) = b.k
get_FEEC_family(::PmLambdaBasis) = :P⁻

function Base.size(b::PmLambdaBasis{D}) where D
  r = get_FEEC_poly_degree(b)
  k = get_FEEC_form_degree(b)
  (binomial(r+k-1,k)*binomial(D+r,D-k), )
end

####################################
# PLambdaBasis nD polynomial bases #
####################################

"""
struct PLambdaBasis{D,V,C,B} <: PolynomialBasis{D,V,Bernstein}

Finite Element Exterior Calculus polynomial basis for the spaces P`ᵣ`Λ`ᴷ` on
`D`-dimensional simplices, but with polynomial forms explicitely transformed
into vectors using the standard equivalence with usual vector calculus (involving
the hodge star operator ⋆ and the flat map ♭).

`V` is `VectorValue{L,T}` where `L=binomial(D,k)`
`C` is binomial(`r`+`k`,`k`)*binomial(`D`+`r`,`D`-`k`), the number of basis polynomials
`B` is the concrete type of the `BernsteinBasisOnSimplex` necessary for the evaluation of the polynomials.

Reference: D. N. Arnold and A. Logg, Periodic Table of the Finite Elements, SIAM News, vol. 47 no. 9, November 2014
"""
struct PLambdaBasis{D,V,C,B} <: PolynomialBasis{D,V,Bernstein}
  r::Int
  k::Int
  scalar_bernstein_basis::B
  Ψ::SVector{C,V}
  _indices::PLambdaIndices

  function PLambdaBasis{D}(::Type{T}, r, k, vertices=nothing;
      diff_geo_calculus_style=false, indices=nothing) where {D,T}

    @check T<:Real "T needs to be <:Real since represents the scalar type, got $T"
    @check k in 0:D "The form order k must be in 0:D, got $k"
    @check r > 0    "The polynomial order r must be positive, got $r"
    if !isnothing(vertices)
      @check length(vertices) == D+1 "$D+1 vertices are required to define a $D-dim simplex, got $(length(vertices))"
      @check eltype(vertices) <: Point{D} "Vertices should be of type <:Point{$D}, got $(eltype(vertices))"
    end

    indices = _generate_PΛ_indices(r,k,D,diff_geo_calculus_style,indices)

    L = binomial(D,k) # Number of components of a basis form
    C = binomial(r+k,k)*binomial(D+r,D-k) # Number of basis polynomials

    if diff_geo_calculus_style
      @notimplemented
    else
      (D>3 && ( 1 < k < D-1)) && @unreachable "Vector calculus proxy of differential form bases not available for D=$D and k=$k"
      V = VectorValue{L,T}
    end

    b = BernsteinBasisOnSimplex{D}(T, r, vertices)
    B = typeof(b)
    Ψ = zero(MVector{C,V})
    _compute_PΛ_basis_form_coefficient!(Ψ,r,k,D,b,vertices,indices)

    if isone(L) && !diff_geo_calculus_style
      V = T
      Ψ = reinterpret(T, Ψ)
    end

    new{D,V,C,B}(r,k,b,Ψ,indices)
  end
end

function PLambdaBasis(::Val{D},::Type{T},r,k,vertices=nothing; diff_geo_calculus_style=false) where {D,T}
  PLambdaBasis{D}(T,r,k,vertices; diff_geo_calculus_style)
end

get_FEEC_poly_degree(b::PLambdaBasis) = b.r
get_FEEC_form_degree(b::PLambdaBasis) = b.k
get_FEEC_family(::PLambdaBasis) = :P

Base.size(::PLambdaBasis{D,V,C}) where {D,V,C} = (C,)


##########################
# Common Implementation  #
##########################

const PΛBases = Union{PmLambdaBasis, PLambdaBasis}

get_order(b::PΛBases) = get_FEEC_poly_degree(b)
get_cart_to_bary_matrix(b::PΛBases) = b.scalar_bernstein_basis.cart_to_bary_matrix
get_bubbles(b::PΛBases) = b._indices.bubbles

function _return_cache(b::PΛBases,x,::Type{G},::Val{N_deriv}) where {G,N_deriv}
  T = eltype(G)
  np = length(x)
  ndof = length(b)
  ndof_bernstein = length(b.scalar_bernstein_basis)

  r = CachedArray(zeros(G,(np,ndof)))
  # Cache for all scalar nD-Bernstein polynomials
  cB = CachedVector(zeros(T,ndof_bernstein))
  if N_deriv > 0
    DB = T
    xi = testitem(x)
    for _ in 1:N_deriv
      DB = gradient_type(DB,xi)
    end
    # Cache for all scalar nD-Bernstein polynomials N_deriv's derivative
    t = (( nothing for _ in 2:N_deriv)..., CachedArray(zeros(DB,(1,ndof_bernstein))))
    s = MArray{Tuple{size(DB)...},T}(undef)
  else
    t = ()
    s = nothing
  end
  (r, s, cB, t...)
end

function _setsize!(b::PΛBases, np, ω, t...)
  ndof = length(b)
  ndof_bernstein = length(b.scalar_bernstein_basis)
  setsize!(ω,(np,ndof))
  setsize!(t[1],(ndof_bernstein,)) # this is cB
  if length(t) > 1
    setsize!(t[end],(1,ndof_bernstein))
  end
end

function _get_static_parameters(b::PΛBases)
  r = get_FEEC_poly_degree(b)
  return Val(r)
end


#################################
# PmLambdaBasis Implementation  #
#################################

function _generate_PmΛ_indices(r,k,D,DG_style,::Nothing)
  identity = objectid( (r,k,D,:P⁻,DG_style) )
  bubbles = PmΛ_bubbles(r,k,D)
  components = basis_forms_components(D,k,DG_style)
  return PLambdaIndices(identity,bubbles,components)
end

function _generate_PmΛ_indices(r,k,D,DG_style,indices::PLambdaIndices)
  @assert objectid( (r,k,D,:P⁻,DG_style) ) == indices.identity
  return indices
end

function _compute_PmΛ_basis_form_coefficient!(m,k::Int,D::Int,b,vertices,indices)
  Vk, VD = Val(k), Val(D)
  _compute_PmΛ_basis_form_coefficient!(m,Vk,VD,b,vertices,indices)
end
function _compute_PmΛ_basis_form_coefficient!(m,Vk,::Val{D},b,vertices,indices) where D
  V = eltype(m)
  M = transpose(b.cart_to_bary_matrix[:,2:end])
  m_J = Mutable(V)(undef)
  @inbounds for (J_id, J) in enumerate(sorted_combinations(Val(D+1),Vk))
    for (I_id, I, I_sgn) in indices.components
      m_J[I_id] = I_sgn*minor(M,I,J,Vk)
    end
    m[J_id] = m_J
  end
  nothing
end

function _compute_PmΛ_basis_form_coefficient!(
  m,Vk::Val{k},::Val{D},b,vertices::Nothing,indices) where {D,k}

  if iszero(k) # so V is scalar, no change of basis
    m .= 1
    return nothing
  end

  V = eltype(m)
  m_J = Mutable(V)(undef)
  @inbounds for (J_id, J) in enumerate(sorted_combinations(Val(D+1),Vk))
    s = Int(isone(J[1]))
    for (I_id, I, I_sgn) in indices.components
      n = count(i-> (J[i]-1)∉I, (1+s):k)
      if iszero(n)
        p = _find_first_val_or_zero(j-> (I[j]+1)∉J, 1, k)
        m_J[I_id] = I_sgn*_minusone_if_odd_else_one(p)
      else
        m_J[I_id] = 0
      end
    end
    m[J_id] = m_J
  end
  nothing
end

# API
function _evaluate_nd!(
  b::PmLambdaBasis{D}, x,
  ω::AbstractMatrix, i, c,
   ::Val{r}) where {D,r}

  V = eltype(ω)
  λ = _cart_to_bary(x, get_cart_to_bary_matrix(b))

  # _evaluate_nd!(::BernsteinBasisOnSimplex) without set_value
  c[1] = 1
  _downwards_de_Casteljau_nD!(c,λ,Val(r-1),Val(D))

  @inbounds for (_, bubble_functions) in get_bubbles(b)
    for (w, _, α_id, J, sub_J_ids) in bubble_functions
      Bα = c[α_id]
      ω_w = zero(V)

      for (l, J_sub_Jl_id) in enumerate(sub_J_ids)
        sgnl = _minusone_if_odd_else_one(l)
        λ_j = λ[J[l]]
        m_J_l = b.m[J_sub_Jl_id]

        ω_w += flipsign(λ_j,sgnl) * m_J_l
      end

      ω[i,w] = Bα * ω_w
    end
  end
end

function _gradient_nd!(
   b::PmLambdaBasis{D}, x,
  ∇ω::AbstractMatrix{G}, i, c,
  ∇B::AbstractMatrix{<:VectorValue{D}},
   s::MVector{D},
   ::Val{r}) where {D,G,r}

  _gradient_nd!(b.scalar_bernstein_basis, x, ∇B, 1, c, nothing, s, Val(r))

  @inbounds for (_, bubble_functions) in get_bubbles(b)
    for (w, α, _, J, sub_J_ids, sup_α_ids) in bubble_functions
      ∇ω_w = zero(G)

      for (l, J_sub_Jl_id) in enumerate(sub_J_ids)
        sgnl = _minusone_if_odd_else_one(l)
        Jl = J[l]
        α_pJl_id = sup_α_ids[Jl]

        ∇Bα_pJl = ∇B[1,α_pJl_id]
        c_α_Jl = (α[Jl]+1) * sgnl / r
        m_J_l = b.m[J_sub_Jl_id]

        ∇ω_w += (c_α_Jl * ∇Bα_pJl) ⊗ m_J_l
      end

      ∇ω[i,w] = ∇ω_w
    end
  end
end

function _hessian_nd!(
   b::PmLambdaBasis{D}, x,
  Hω::AbstractMatrix{G}, i, c,
    ::Nothing,
  HB::AbstractMatrix{<:TensorValue{D,D}},
   s::MMatrix{D,D},
   ::Val{r}) where {D,G,r}

  _hessian_nd!(b.scalar_bernstein_basis, x, HB, 1, c, nothing, nothing, s, Val(r))

  @inbounds for (_, bubble_functions) in get_bubbles(b)
    for (w, α, _, J, sub_J_ids, sup_α_ids) in bubble_functions
      Hω_w = zero(G)

      for (l, J_sub_Jl_id) in enumerate(sub_J_ids)
        sgnl = _minusone_if_odd_else_one(l)
        Jl = J[l]
        α_pJl_id = sup_α_ids[Jl]

        HB_αJl = HB[1,α_pJl_id]
        c_αJl = (α[Jl]+1) * sgnl / r
        m_Jl = b.m[J_sub_Jl_id]

        Hω_w += (c_αJl * HB_αJl) ⊗ m_Jl
      end

      Hω[i,w] = Hω_w
    end
  end
end


#################################
# PLambdaBasis Implementation  #
#################################

function _generate_PΛ_indices(r,k,D,DG_style,::Nothing)
  identity = objectid( (r,k,D,:P,DG_style) )
  bubbles = PΛ_bubbles(r,k,D)
  components = basis_forms_components(D,k,DG_style)
  return PLambdaIndices(identity,bubbles,components)
end

function _generate_PΛ_indices(r,k,D,DG_style,indices::PLambdaIndices)
  @assert objectid( (r,k,D,:P,DG_style) ) == indices.identity
  return indices
end

function _compute_PΛ_basis_form_coefficient!(Ψ,r::Int,k::Int,D::Int,b,vertices,indices)
  Vr, Vk, VD = Val(r), Val(k), Val(D)
  _compute_PΛ_basis_form_coefficient!(Ψ,Vr,Vk,VD,b,vertices,indices)
end
function _compute_PΛ_basis_form_coefficient!(Ψ,Vr,Vk,::Val{D},b,vertices,indices) where D
  N = D+1
  V = eltype(Ψ)
  T = eltype(V)
  α_prec = ntuple(_->-1, N)
  φ_αF = MMatrix{D,N,T}(undef)
  Ψw = Mutable(V)(undef)
  @inbounds for (F, bubble_functions) in indices.bubbles
    for (w, α, _, J) in bubble_functions
      if α ≠ α_prec
        update_φ_αF!(φ_αF,b,α,F,Vr)
        α_prec = α
      end

      for (I_id, I, I_sgn) in indices.components
        Ψw[I_id] = I_sgn * minor(φ_αF,I,J,Vk)
      end
      Ψ[w] = Ψw
    end
  end
  nothing
end

@inline function update_φ_αF!(φ_αF,b,α,F,::Val{r}) where r
  M = b.cart_to_bary_matrix
  @inbounds for ci in CartesianIndices(φ_αF)
    i, j = ci[1], ci[2]
    mF = sum(M[Fl,i+1] for Fl in F; init=0)
    φ_αF[ci] = M[j,i+1] - α[j]*mF/r
  end
end

function _compute_PΛ_basis_form_coefficient!(
  Ψ,Vr,Vk,::Val,b,vertices::Nothing,indices)

  V = eltype(Ψ)
  T = eltype(V)
  Ψw = Mutable(V)(undef)
  @inbounds for (F, bubble_functions) in indices.bubbles
    for (w, α, _, J) in bubble_functions
      for (I_id, I, I_sgn) in indices.components
        Ψw[I_id] = I_sgn * _hat_Ψ(Vr,Vk,α,F,I,J,T)
      end
      Ψ[w] = Ψw
    end
  end
  nothing
end

"""
    _hat_Ψ(::Val{r},::Val{k},α,F,I,J,T)::T

PLambdaBasis.Ψ matrix elements in the reference simplex, T is the scalar return type

This is actually not faster than computing the matrices and the minors
explicitely like when vertices are given, but lets keep it here in case we want
to compute these at compile time one day.
"""
function _hat_Ψ(Vr::Val{r},Vk::Val{k},α,F,I,J,::Type{T})::T where {T,r,k}
  @check sum(α) == r
  @check length(I) == length(J) == k

  iszero(k) && return one(T) # 0 forms

  @inbounds begin

    s = Int(isone(J[1]))
    n = count(i-> (J[i]-1)∉I, (1+s):k)

    n > 1 && return 0. # rank M_IJ inferior to 2

    p = _find_first_val_or_zero(j-> (I[j]+1)∉J, 1, k)

    if isone(n)        # rank M_IJ is 1
      m = _find_first_val_or_zero(i-> (J[i]-1)∉I, (s+1), k)
      u_p, v_m = _u(p,F,I), _v(m,α,J,Vr)
      sgn = isodd(m+p) ? -1 : 1
      iszero(s) && return sgn*u_p*v_m

      q = _find_first_val_or_zero(j-> (I[j]+s)∉J, (p+1), k)
      u_q = _u(q,F,I)
      sgn *= isodd(q) ? -1 : 1
      return sgn * v_m * (u_q - u_p)
    end

    u, v = _u(F,I,Vk), _v(α,J,Vr,Vk)
    if iszero(s)
      return 1 + sum( u .* v )
    else
      Ψ_IJ = one(T)
      sum_v = v[1]
      for l in 1:p-1
        vlp = v[l+1]
        sum_v += vlp
        Ψ_IJ += vlp*u[l]
      end
      for l in (p+1):k
        vl = v[l]
        sum_v += vl
        Ψ_IJ += vl*u[l]
      end
      sgn = _minusone_if_odd_else_one(p)
      return sgn * (Ψ_IJ - u[p]*sum_v)
    end

  end
  @unreachable
end

@propagate_inbounds _u(i::Int,F,I)   = Int(isone(F[1])) - Int(I[i]+1 in F)
@propagate_inbounds _u(F::Combination,I,Vk) = ntuple(i->_u(i,F,I), Vk)
@propagate_inbounds _v(j::Int,α,J,::Val{r}) where r = α[J[j]]/r
@propagate_inbounds _v(α::BernsteinTerm,J,Vr,Vk) = ntuple(j->_v(j,α,J,Vr), Vk)


# API

function _evaluate_nd!(
  b::PLambdaBasis{D}, x,
  ω::AbstractMatrix, i, c,
   ::Val{r}) where {D,r}

  λ = _cart_to_bary(x, get_cart_to_bary_matrix(b))

  # _evaluate_nd!(::BernsteinBasisOnSimplex) without set_value
  c[1] = 1
  _downwards_de_Casteljau_nD!(c,λ,Val(r),Val(D))

  @inbounds for (_, bubble_functions) in get_bubbles(b)
    for (w, _, α_id) in bubble_functions
      Ψw = b.Ψ[w]
      Bα = c[α_id]
      ω[i,w] = Bα * Ψw
    end
  end
end

function _gradient_nd!(
   b::PLambdaBasis{D}, x,
  ∇ω::AbstractMatrix, i, c,
  ∇B::AbstractMatrix{<:VectorValue{D}},
   s::MVector{D},
   ::Val{r}) where {D,r}

  _gradient_nd!(b.scalar_bernstein_basis, x, ∇B, 1, c, nothing, s, Val(r))

  @inbounds for (_, bubble_functions) in get_bubbles(b)
    for (w, _, α_id) in bubble_functions
      Ψw = b.Ψ[w]
      ∇Bα = ∇B[1,α_id]
      ∇ω[i,w] = ∇Bα ⊗ Ψw
    end
  end
end

function _hessian_nd!(
   b::PLambdaBasis{D}, x,
  Hω::AbstractMatrix, i, c,
    ::Nothing,
  HB::AbstractMatrix{<:TensorValue{D,D}},
   s::MMatrix{D,D},
   ::Val{r}) where {D,r}

  _hessian_nd!(b.scalar_bernstein_basis, x, HB, 1, c, nothing, nothing, s, Val(r))

  @inbounds for (_, bubble_functions) in get_bubbles(b)
    for (w, _, α_id) in bubble_functions
      Ψw = b.Ψ[w]
      HBα = HB[1,α_id]
      Hω[i,w] = HBα ⊗ Ψw
    end
  end
end


##########################
# PLambda bases helpers  #
##########################

@propagate_inbounds function _minusone_if_odd_else_one(i)
  isodd(i) ? -1 : 1
end
@propagate_inbounds function _find_first_val_or_zero(pred, start, stop)
  r = findfirst(pred,start:stop)
  return isnothing(r) ? 0 : r+start-1
end

function _PmΛ_F_bubble_functions(r,k,D,F,i)
  N = D + 1
  ids = BubbleFunction[]
  for α in bernstein_terms(r-1,N)
    α = Int[α...]
    sup_α_ids = _sup_multi_indices(α)
    for J in sorted_combinations(N,k+1)
      J = Int[J...]
      sub_J_ids = _sub_combinations_ids(J)
      j = minimum(J)-1
      if issetequal(support(α) ∪ J, F) && all(α[1:j] .== 0)
        i += 1
        α_id = bernstein_term_id(α)
        push!(ids, (i, α, α_id, J, sub_J_ids, sup_α_ids))
      end
    end
  end
  ids
end

function PmΛ_bubbles(r,k,D)
  i=0
  F_bubble_functions = Bubble[]
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      F = Int[F...]
      bubble_functions = _PmΛ_F_bubble_functions(r,k,D,F,i)
      push!(F_bubble_functions, (F, bubble_functions))
      i += length(bubble_functions)
    end
  end
  @check i == binomial(r+k-1,k)*binomial(D+r,D-k)
  F_bubble_functions
end

function _PΛ_F_bubble_functions(r,k,D,F,i)
  N = D + 1
  ids = BubbleFunction[]
  empty_vec = Int[]
  for α in bernstein_terms(r,D)
    α = Int[α...]
    for J in sorted_combinations(N,k)
      J = Int[J...]
      j = minimum(setdiff(F,J), init=N+1)
      j = j==(N+1) ? 0 : j-1
      if issetequal(support(α) ∪ J, F) && all(α[1:j] .== 0)
        i += 1
        α_id = bernstein_term_id(α)
        push!(ids, (i, α, α_id, J, empty_vec, empty_vec))
      end
    end
  end
  ids
end

function PΛ_bubbles(r,k,D)
  i=0
  F_bubble_functions = Bubble[]
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      F = Int[F...]
      bubble_functions = _PΛ_F_bubble_functions(r,k,D,F,i)
      push!(F_bubble_functions, (F, bubble_functions))
      i += length(bubble_functions)
    end
  end
  @check i == binomial(r+k,k)*binomial(D+r,D-k)
  F_bubble_functions
end

function basis_forms_components(D,k,DG_style)
  components = Vector{Component}(undef, binomial(D,k))
  for (I_id, I) in enumerate(sorted_combinations(D,k))
    I = Int[I...]
    if DG_style
      components[I_id] = (I_id, I, 1)
    else
      Icomp = complement(I, D)
      Istar_id = combination_index(Icomp)
      Icomp_sgn = combination_sign(I)
      components[I_id] = (Istar_id, I, Icomp_sgn)
    end
  end
  components
end

