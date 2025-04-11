#####################################
# P⁻LambdaBasis nD polynomial bases #
#####################################

"""
struct PmLambdaBasis{D,V,L,B} <: PolynomialBasis{D,V,Bernstein}

Finite Element Exterior Calculus polynomial basis for the spaces `D`-dimensional
simplices, P`r`⁻Λ`ᴷ`(△`ᴰ`), but with polynomial forms explicitely transformed
into vectors using the standard equivalence with usual vector calculus (given
by the flat map ♭).

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
  indices::PLambdaIndices

  function PmLambdaBasis{D}(::Type{T}, r, k, vertices=nothing;
      diff_geo_calculus_style=false,
      indices=nothing) where {D,T}

    @check T<:Real "T needs to be <:Real since represents the scalar type, got $T"
    @check k in 0:D "The form order k must be in 0:D, got $k"
    @check r > 0    "The polynomial order r must be positive, got $r"
    if !isnothing(vertices)
      @check length(vertices) == D+1 "$D+1 vertices are required to define a $D-dim simplex, got $(length(vertices))"
      @check eltype(vertices) <: Point{D} "Vertices should be of type <:Point{$D}, got $(eltype(vertices))"
    end

    pmΛ_indices = _generate_PmΛ_indices(r,k,D,:P⁻,indices)

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
    _compute_PmΛ_basis_form_coefficient!(m,k,D,b,vertices,pmΛ_indices)

    if isone(L) && !diff_geo_calculus_style
      V = T
      m = reinterpret(T, m)
    end

    new{D,V,LN,B}(r,k,b,m,pΛ_indices)
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
  binomial(r+k,k)*binomial(D+r,D-k)
end

####################################
# PLambdaBasis nD polynomial bases #
####################################

"""
struct PLambdaBasis{D,V,C,B} <: PolynomialBasis{D,V,Bernstein}

Finite Element Exterior Calculus polynomial basis for the spaces `D`-dimensional
simplices, P`r`Λ`ᴷ`(△`ᴰ`), but with polynomial forms explicitely transformed
into vectors using the standard equivalence with usual vector calculus (given
by the flat map ♭).

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
  indices::PLambdaIndices

  function PLambdaBasis{D}(::Type{T}, r, k, vertices=nothing;
      diff_geo_calculus_style=false,
      pΛ_indices=nothing) where {D,T}

    @check T<:Real "T needs to be <:Real since represents the scalar type, got $T"
    @check k in 0:D "The form order k must be in 0:D, got $k"
    @check r > 0    "The polynomial order r must be positive, got $r"
    if !isnothing(vertices)
      @check length(vertices) == D+1 "$D+1 vertices are required to define a $D-dim simplex, got $(length(vertices))"
      @check eltype(vertices) <: Point{D} "Vertices should be of type <:Point{$D}, got $(eltype(vertices))"
    end

    pΛ_indices = _generate_PΛ_indices(r,k,D,pΛ_indices)

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
    _compute_PΛ_basis_form_coefficient!(Ψ,r,k,D,b,vertices,pΛ_indices)

    if isone(L) && !diff_geo_calculus_style
      V = T
      Ψ = reinterpret(T, Ψ)
    end

    new{D,V,C,B}(r,k,b,Ψ,pΛ_indices)
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
get_bubbles(b::PΛBases) = b.indices.bubbles

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
  setsize!(t[1],(ndof_bernstein,))
  if length(t) > 1
    setsize!(t[end],(1,ndof_bernstein))
  end
end

function _get_static_parameters(b::PΛBases)
  r = get_FEEC_poly_degree(b)
  k = get_FEEC_form_degree(b)
  return (Val(r), Val(k))
end


#################################
# PmLambdaBasis Implementation  #
#################################

function _generate_PmΛ_indices(r,k,D,P,::Nothing)
  identity = objectid( (r,k,D,P) )
  bubbles = PmΛ_bubbles(r,k,D)
  # if _DG_calculus_style(V) end # false ATM
  components = collect(enumerate(sorted_combinations(D,k))) #TODO dpd on style
  b = Vector{Tuple{Int, Int, Combination}}[]
  return PLambdaIndices(identity,bubbles,components,b)
end

function _generate_PmΛ_indices(r,k,D,P,pΛ_indices)
  @assert objectid( (r,k,D,:P⁻) ) == pΛ_indices.identity
  return pΛ_indices
end

function _compute_PmΛ_basis_form_coefficient!(m,k::Int,D::Int,b,vertices,pmΛ_indices)
  Vk, VD = Val(k), Val(D)
  _compute_PmΛ_basis_form_coefficient!(m,Vk,VD,b,vertices,pmΛ_indices)
end
function _compute_PmΛ_basis_form_coefficient!(m,Vk,VD::Val{D},b,vertices,pmΛ_indices) where D
  V = eltype(m)
  M = transpose(b.cart_to_bary_matrix[:,2:end])
  m_J = Mutable(V)(undef)
  @inbounds for (J_id,J) in enumerate(sorted_combinations(Val(D+1),Vk))
    for (I_id,I) in pmΛ_indices.components
      m_J[I_id] = minor(M,I,J,Vk)
    end
    m[J_id] = m_J
  end
  nothing
end

function _compute_PmΛ_basis_form_coefficient!(
  m,Vk::Val{k},VD::Val{D},b,vertices::Nothing,pmΛ_indices) where {D,k}

  if iszero(k) # so V is scalar, no change of basis
    m .= 1
    return nothing
  end

  V = eltype(m)
  m_J = Mutable(V)(undef)
  @inbounds for (J_id,J) in enumerate(sorted_combinations(Val(D+1),Vk))
    s = Int(isone(J[1]))
    for (I_id,I) in pmΛ_indices.components
      n = count(i-> (J[i]-1)∉I, (1+s):k)
      if iszero(n)
        p = _find_first_val_or_zero(j-> (I[j]+1)∉J, 1, k)
        m_J[I_id] = _minusone_if_odd_else_one(p)
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
  b::PmLambdaBasis{D,V}, x,
  ω::AbstractMatrix{V}, i, c,
   ::Tuple{Val{r},Val{k}}) where {D,V,r,k}

  λ = _cart_to_bary(x, get_cart_to_bary_matrix(b))

  # _evaluate_nd!(::BernsteinBasisOnSimplex) without set_value
  c[1] = 1
  _downwards_de_Casteljau_nD!(c,λ,Val(r-1),Val(D))

  sub_comb_ids = MVector{k+1,Tuple{Int,Int}}(undef) # J is of length k+1

  @inbounds for (_, _, dF_bubbles) in get_bubbles(b)#PmΛ_bubbles(Val(r),Val(k),Val(D))
    for (w, _, α_id, J) in dF_bubbles
      Bα = c[α_id]
      ω_w = zero(V)

      _sub_combinations_ids!(sub_comb_ids, J)
      for (l,J_sub_l_id) in sub_comb_ids
        sgn = _minusone_if_odd_else_one(l)
        λ_j = λ[J[l]]
        m_Jl = b.m[J_sub_l_id]

        ω_w += sgn * λ_j * m_Jl
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
    ::Tuple{Val{r},Val{k}}) where {D,G,r,k}

  _gradient_nd!(b.scalar_bernstein_basis, x, ∇B, 1, c, nothing, s, Val(r))

  sub_comb_ids = MVector{k+1,Tuple{Int,Int}}(undef) # J is of length k+1

  @inbounds for (_, _, dF_bubbles) in get_bubbles(b)#PmΛ_bubbles(Val(r),Val(k),Val(D))
    for (w, α, _, J) in dF_bubbles
      ∇ω_w = zero(G)

      _sub_combinations_ids!(sub_comb_ids, J)
      for (l,J_sub_l_id) in sub_comb_ids
        j = J[l]
        α_pj = α .+ ntuple(i->Int(i==j), Val(D+2))
        α_pj_id = bernstein_term_id(α_pj)

        ∇Bα_pj = ∇B[1,α_pj_id]
        s_α_j = α_pj[j] * _minusone_if_odd_else_one(l)
        m_Jl = b.m[J_sub_l_id]

        ∇ω_w += (s_α_j .* ∇Bα_pj) ⊗ m_Jl
      end

      ∇ω[i,w] = ∇ω_w / r
    end
  end
end

function _hessian_nd!(
   b::PmLambdaBasis{D}, x,
  Hω::AbstractMatrix{G}, i, c,
    ::Nothing,
  HB::AbstractMatrix{<:TensorValue{D,D}},
   s::MMatrix{D,D},
    ::Tuple{Val{r},Val{k}}) where {D,G,r,k}

  _hessian_nd!(b.scalar_bernstein_basis, x, HB, 1, c, nothing, nothing, s, Val(r))

  sub_comb_ids = MVector{k+1,Tuple{Int,Int}}(undef) # J is of length k+1

  @inbounds for (_, _, dF_bubbles) in get_bubbles(b)#PmΛ_bubbles(Val(r),Val(k),Val(D))
    for (w, α, _, J) in dF_bubbles
      Hω_w = zero(G)

      _sub_combinations_ids!(sub_comb_ids, J)
      for (l,J_sub_l_id) in sub_comb_ids
        j = J[l]
        α_pj = α .+ ntuple(i->Int(i==j), Val(D+2))
        α_pj_id = bernstein_term_id(α_pj)

        HBα_pj = HB[1,α_pj_id]
        s_α_j = α_pj[j] * _minusone_if_odd_else_one(l)
        m_Jl = b.m[J_sub_l_id]

        Hω_w +=  (s_α_j .* HBα_pj) ⊗ m_Jl
      end

      Hω[i,w] = Hω_w / r
    end
  end
end


#################################
# PLambdaBasis Implementation  #
#################################

function _generate_PΛ_indices(r,k,D,::Nothing)
  identity = objectid( (r,k,D,:P) )
  bubbles = PΛ_bubbles(r,k,D)
  # if _DG_calculus_style(V) end # false ATM
  components = enumerate(sorted_combinations(D,k))#TODO dpd on style
  b = Vector{Tuple{Int, Int, Combination}}[]
  return PLambdaIndices(identity,bubbles,components,b)
end

function _generate_PΛ_indices(r,k,D,pΛ_indices)
  @assert objectid( (r,k,D,:P) ) == pΛ_indices.identity
  return pΛ_indices
end

function _compute_PΛ_basis_form_coefficient!(Ψ,r::Int,k::Int,D::Int,b,vertices,pΛ_indices)
  Vr, Vk, VD = Val(r), Val(k), Val(D)
  _compute_PΛ_basis_form_coefficient!(Ψ,Vr,Vk,VD,b,vertices,pΛ_indices)
end
function _compute_PΛ_basis_form_coefficient!(Ψ,Vr,Vk,VD::Val{D},b,vertices,pΛ_indices) where D
  N = D+1
  V = eltype(Ψ)
  T = eltype(V)
  α_prec = ntuple(_->-1, N)
  φ_αF = MMatrix{D,N,T}(undef)
  Ψw = Mutable(V)(undef)
  @inbounds for (_, F, dF_bubbles) in pΛ_indices.bubbles
    for (w, α, _, J) in dF_bubbles
      if α ≠ α_prec
        update_φ_αF!(φ_αF,b,α,F,Vr)
        α_prec = α
      end

      for (sgnIcomp, Istar_id, I) in pΛ_indices.components
        Ψw[Istar_id] = sgnIcomp * minor(φ_αF,I,J,Vk)
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
  Ψ,Vr,Vk,VD::Val{D},b,vertices::Nothing,pΛ_indices) where D

  V = eltype(Ψ)
  T = eltype(V)
  Ψw = Mutable(V)(undef)
  @inbounds for (_, F, dF_bubbles) in pΛ_indices.bubbles
    for (w, α, _, J) in dF_bubbles
      for (Istar_id, I, sgnIcomp) in pΛ_indices.components
        Ψw[Istar_id] = sgnIcomp * _hat_Ψ(Vr,Vk,α,F,I,J,T)
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
  b::PLambdaBasis{D,V}, x,
  ω::AbstractMatrix{V}, i, c,
   ::Tuple{Val{r},Val{k}}) where {D,V,r,k}

  λ = _cart_to_bary(x, get_cart_to_bary_matrix(b))

  # _evaluate_nd!(::BernsteinBasisOnSimplex) without set_value
  c[1] = 1
  _downwards_de_Casteljau_nD!(c,λ,Val(r),Val(D))

  @inbounds for (_, _, dF_bubbles) in get_bubbles(b)#PΛ_bubbles(Val(r),Val(k),Val(D))
    for (w, _, α_id, _) in dF_bubbles
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
    ::Tuple{Val{r},Val{k}}) where {D,r,k}

  _gradient_nd!(b.scalar_bernstein_basis, x, ∇B, 1, c, nothing, s, Val(r))

  @inbounds for (_, _, dF_bubbles) in get_bubbles(b)#PΛ_bubbles(Val(r),Val(k),Val(D))
    for (w, _, α_id, _) in dF_bubbles
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
    ::Tuple{Val{r},Val{k}}) where {D,r,k}

  _hessian_nd!(b.scalar_bernstein_basis, x, HB, 1, c, nothing, nothing, s, Val(r))

  @inbounds for (_, _, dF_bubbles) in get_bubbles(b)#PΛ_bubbles(Val(r),Val(k),Val(D))
    for (w, _, α_id, _) in dF_bubbles
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

function _PmΛ_F_bubbles(r,k,D,F,i)
  N = D + 1
  ids = Tuple{Int, BernsteinTerm, Int, Combination}[]
  for α in bernstein_terms(r-1,N)
    for J in sorted_combinations(N,k+1)
      J = Int[J...]
      j = minimum(J)-1
      if issetequal(supp(α) ∪ J, F) && all(α[1:j] .== 0)
        i += 1
        α_id = bernstein_term_id(α)
        push!(ids, (i, copy(α), α_id, J))
      end
    end
  end
  ids
end

PmΛ_bubbles(r,k,D) = PmΛ_bubbles(Val(r),Val(k),Val(D))
@generated function PmΛ_bubbles(::Val{r},::Val{k},::Val{D}) where {r,k,D}
  i=0
  d_F_bubbles = Tuple{Int, Combination, Vector{Tuple{Int, BernsteinTerm, Int, Combination}}}[]
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      F = Int[F...]
      dF_bubbles = _PmΛ_F_bubbles(r,k,D,F,i)
      push!(d_F_bubbles, (d, F, dF_bubbles))
      i += length(dF_bubbles)
    end
  end
  @check i == binomial(r+k-1,k)*binomial(D+r,D-k)
  :( $(d_F_bubbles) )
end

function _PΛ_F_bubbles(r,k,D,F,i)
  N = D + 1
  ids = Tuple{Int, BernsteinTerm, Int, Combination}[]
  for α in bernstein_terms(r,D)
    for J in sorted_combinations(N,k)
      J = Int[J...]
      j = minimum(setdiff(F,J), init=N+1)
      j = j==(N+1) ? 0 : j-1
      if issetequal(supp(α) ∪ J, F) && all(α[1:j] .== 0)
        i += 1
        α_id = bernstein_term_id(α)
        push!(ids, (i, copy(α), α_id, J))
      end
    end
  end
  ids
end

PΛ_bubbles(r,k,D) = PΛ_bubbles(Val(r),Val(k),Val(D))
@generated function PΛ_bubbles(::Val{r},::Val{k},::Val{D}) where {r,k,D}
  i=0
  d_F_bubbles = Tuple{Int, Combination, Vector{Tuple{Int, BernsteinTerm, Int, Combination}}}[]
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      F = Int[F...]
      dF_bubbles = _PΛ_F_bubbles(r,k,D,F,i)
      push!(d_F_bubbles, (d, F, dF_bubbles))
      i += length(dF_bubbles)
    end
  end
  @check i == binomial(r+k,k)*binomial(D+r,D-k)
  :( $(d_F_bubbles) )
end

#_all_combis_and_sub_combi_ids(::Val,::Val{0}) = Tuple{}[]
#@generated function _all_combis_and_sub_combi_ids(::Val{D},::Val{k}) where {D,k}
#  @check k>0
#
#  res = Tuple{Int, Combination, Vector{Tuple{Int, Int}}}[]
#  sub_comb_ids = MVector{k+1,Tuple{Int,Int}}(undef) # J is of length k+1
#
#  for (I_id,I) in enumerate(sorted_combinations(Val(D),Val(k)))
#    _sub_combinations_ids!(sub_comb_ids, I)
#    push!(res, (I_id, I, sub_comb_ids))
#  end
#  :( $(res) )
#end

@generated function hodged_basis_forms(::Val{D},::Val{k}) where {D,k}
  res = Tuple{Int,Int,Combination}[]
  for I in sorted_combinations(Val(D),Val(k))
    I = Int[I...]
    Icomp = complement(I, Val(D))
    Istar_id = combination_index(Icomp)
    sgnIcomp = combination_sign(I)
    # actually ⋆(I) = Istar = sgnIcomp*Icomp
    push!(res, (sgnIcomp, Istar_id, I) )
  end
  :( $(res) )
end
