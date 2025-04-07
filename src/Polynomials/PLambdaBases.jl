#################################
# Tensorial nD polynomial bases #
#################################

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

  function PLambdaBasis{D}(::Type{T}, r, k, vertices=nothing; diff_geo_calculus_style=false) where {D,T}
    @check T<:Real "T needs to be <:Real since represents the scalar type"
    @check k in 0:D "The form order k must be in 0:D"
    @check r > 0    "The polynomial order r must be positive"
    if !isnothing(vertices)
      @check length(vertices) == D+1 "$D+1 vertices are required to define a $D-dim simplex, got $(length(vertices))"
      @check eltype(vertices) <: Point{D} "Vertices should be of type Point{$D}, got $(eltype(vertices))"
    end

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
    _compute_PΛ_basis_form_coefficient!(Ψ,r,k,D,b,vertices)

    if isone(L) && !diff_geo_calculus_style
      V = T
      Ψ = reinterpret(T, Ψ)
    end

    new{D,V,C,B}(r,k,b,Ψ)
  end
end

function PLambdaBasis(::Val{D},::Type{T},r,k,vertices=nothing; diff_geo_calculus_style=false) where {D,T}
  PLambdaBasis{D}(T,r,k,vertices; diff_geo_calculus_style)
end

"""
    _DG_calculus_style(::V) = false

Temporary API to signal that the coefficient of `V`-valued polynomial forms
should be transformed into classic vector calulus components using Hodge star
operator.
"""
_DG_calculus_style(::V) where V = false

get_FEEC_poly_degree(b::PLambdaBasis) = b.r
get_FEEC_form_degree(b::PLambdaBasis) = b.k
get_FEEC_family(::PLambdaBasis) = :P

Base.size(::PLambdaBasis{D,V,C}) where {D,V,C} = (C,)
get_order(b::PLambdaBasis) = get_FEEC_poly_degree(b)

get_cart_to_bary_matrix(b::PLambdaBasis) = b.scalar_bernstein_basis.cart_to_bary_matrix
#get_dimension(::PLambdaBasis{D}) where D = D


###################
# Implementation  #
###################

function _compute_PΛ_basis_form_coefficient!(Ψ,r::Int,k::Int,D::Int,b,vertices)
  Vr, Vk, VD = Val(r), Val(k), Val(D)
  _compute_PΛ_basis_form_coefficient!(Ψ,Vr,Vk,VD,b,vertices)
end
function _compute_PΛ_basis_form_coefficient!(Ψ,Vr,Vk,VD::Val{D},b,vertices) where D
  N = D+1
  V = eltype(Ψ)
  T = eltype(V)
  α_prec = ntuple(_->-1, N)
  φ_αF = MMatrix{D,N,T}(undef)
  Ψw = Mutable(V)(undef)
  @inbounds for (_, F, dF_bubbles) in PΛ_bubble_indices(Vr,Vk,VD)
    for (w, α, _, J) in dF_bubbles
      if α ≠ α_prec
        update_φ_αF!(φ_αF,b,α,F,Vr)
        α_prec = α
      end
      if _DG_calculus_style(V) # false ATM
        for (I_id,I) in enumerate(sorted_combinations(VD,Vk))
          Ψw[I_id] = @inline minor(φ_αF,I,J)
        end
      else
        for (sgnIcomp, Istar_id, I) in hodged_basis_forms(VD,Vk)
          Ψw[Istar_id] = sgnIcomp * @inline minor(φ_αF,I,J)
        end
      end
      Ψ[w] = Ψw
    end
  end
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
  Ψ,Vr,Vk,VD::Val{D},b,vertices::Nothing) where D

  V = eltype(Ψ)
  T = eltype(V)
  Ψw = Mutable(V)(undef)
  @inbounds for (_, F, dF_bubbles) in PΛ_bubble_indices(Vr,Vk,VD)
    for (w, α, _, J) in dF_bubbles
      if _DG_calculus_style(V) # false ATM
        for (I_id,I) in enumerate(sorted_combinations(VD,Vk))
          Ψw[I_id] = _hat_Ψ(Vr,α,F,I,J,T)
        end
      else
        for (sgnIcomp, Istar_id, I) in hodged_basis_forms(VD,Vk)
          Ψw[Istar_id] = sgnIcomp * _hat_Ψ(Vr,α,F,I,J,T)
        end
      end
      Ψ[w] = Ψw
    end
  end
end

"""
    _hat_Ψ(::Val{r},α,F,I,J,T)::T

PLambdaBasis.Ψ matrix elements in the reference simplex, T is the scalar return type
"""
function _hat_Ψ(Vr,α,F,I,J,::Type{T})::T where T
  @check Val(sum(α)) == Vr
  @check length(I) == length(J) # thanks to dispatch on zero length J below

  iszero(length(J)) && return one(T) # 0 forms
  Vk = Val(length(J))

  @inbounds begin

    s = Int(isone(J[1]))
    n = count(i-> (J[i]-1)∉I, (s+1):length(J))

    n > 1 && return 0. # rank M_IJ inferior to 2

    p = _find_first_val_or_zero(j-> (I[j]+1)∉J, 1, length(J))

    if isone(n)        # rank M_IJ is 1
      m = _find_first_val_or_zero(i-> (J[i]-1)∉I, (s+1), length(J))
      u_m, v_p = _u(m,F,I), _v(p,α,J,Vr)
      iszero(s) && return (-1)^(m+p)*u_m*v_p

      q = _find_first_val_or_zero(j-> (I[j]+s)∉J, (p+1), length(J))
      @check !iszero(q)
      v_q = _v(q,α,J,Vr)
      return (-1)^(m+p+q) * u_m * (v_q - v_p)
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
      for l in (p+1):length(J)
        vl = v[l]
        sum_v += vl
        Ψ_IJ += vl*u[l]
      end
      return (-1)^p * (Ψ_IJ - u[p]*sum_v)
    end

  end
  @unreachable
end

@propagate_inbounds function _find_first_val_or_zero(pred, start, stop)
  r = findfirst(pred,start:stop)
  return isnothing(r) ? 0 : r+start-1
end

@propagate_inbounds _u(i::Int,F,I)   = Int(isone(F[1])) - Int(I[i]+1 in F)
@propagate_inbounds _u(F::Combination,I,Vk) = ntuple(i->_u(i,F,I), Vk)
@propagate_inbounds _v(j::Int,α,J,::Val{r}) where r = α[J[j]]/r
@propagate_inbounds _v(α::Tuple,J,Vr,Vk) = ntuple(j->_v(j,α,J,Vr), Vk)


# API

function _return_cache(
  b::PLambdaBasis{D},x,::Type{G},::Val{N_deriv}) where {D,G,N_deriv}

  @assert D == length(eltype(x)) "Incorrect number of point components"
  T = eltype(G)
  r = get_order(b)
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

function _setsize!(b::PLambdaBasis{D}, np, ω, t...) where D
  r = get_order(b)
  ndof = length(b)
  ndof_bernstein = length(b.scalar_bernstein_basis)
  setsize!(ω,(np,ndof))
  setsize!(t[1],(ndof_bernstein,))
  if length(t) > 1
    setsize!(t[end],(1,ndof_bernstein))
  end
end

function _get_static_parameters(b::PLambdaBasis)
  r = get_FEEC_poly_degree(b)
  k = get_FEEC_form_degree(b)
  return (Val(r), Val(k))
end

function _evaluate_nd!(
  b::PLambdaBasis{D,V}, x,
  ω::AbstractMatrix{V}, i,
  c::AbstractVector{T},
  ::Tuple{Val{r},Val{k}}) where {D,V,T,r,k}

  λ = _cart_to_bary(x, get_cart_to_bary_matrix(b))

  # _evaluate_nd!(::BernsteinBasisOnSimplex) without set_value
  c[1] = one(T)
  _downwards_de_Casteljau_nD!(c,λ,Val(r),Val(D))

  @inbounds for (_, _, dF_bubbles) in PΛ_bubble_indices(Val(r),Val(k),Val(D))
    for (w, _, α_id, _) in dF_bubbles
      Ψw = b.Ψ[w]
      Bα = c[α_id]
      ω[i,w] = Bα * Ψw
    end
  end
end

function _gradient_nd!(
  b::PLambdaBasis{D}, x,
  ∇ω::AbstractMatrix{G}, i,
   c::AbstractVector{T},
  ∇B::AbstractMatrix{<:VectorValue{D,T}},
   s::MVector{D,T},
  ::Tuple{Val{r},Val{k}}) where {D,G,T,r,k}

  _gradient_nd!(b.scalar_bernstein_basis, x, ∇B, 1, c, nothing, s, Val(r))

  @inbounds for (_, _, dF_bubbles) in PΛ_bubble_indices(Val(r),Val(k),Val(D))
    for (w, _, α_id, _) in dF_bubbles
      Ψw = b.Ψ[w]
      ∇Bα = ∇B[1,α_id]
      ∇ω[i,w] = ∇Bα ⊗ Ψw
    end
  end
end

function _hessian_nd!(
  b::PLambdaBasis{D}, x,
  Hω::AbstractMatrix{G}, i,
   c::AbstractVector{T},
    ::Nothing,
  HB::AbstractMatrix{<:TensorValue{D,D,T}},
   s::MMatrix{D,D,T},
  ::Tuple{Val{r},Val{k}}) where {D,G,T,r,k}

  _hessian_nd!(b.scalar_bernstein_basis, x, HB, 1, c, nothing, nothing, s, Val(r))

  @inbounds for (_, _, dF_bubbles) in PΛ_bubble_indices(Val(r),Val(k),Val(D))
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

function _P⁻Λ_F_bubble_indices(r,k,D,F,i)
  N = D + 1
  ids = Tuple{Int64, BernsteinTerm, Int, Combination{N}}[]
  for α in bernstein_terms(r-1,N)
    for J in sorted_combinations(N,k+1)
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

function _PΛ_F_bubble_indices(r,k,D,F,i)
  N = D + 1
  ids = Tuple{Int64, BernsteinTerm, Int, Combination{N}}[]
  for α in bernstein_terms(r,D)
    for J in sorted_combinations(N,k)
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

PΛ_bubble_indices(r,k,D) = PΛ_bubble_indices(Val(r),Val(k),Val(D))
@generated function PΛ_bubble_indices(::Val{r},::Val{k},::Val{D}) where {r,k,D}
  i=0
  N = D+1
  d_F_bubbles = Tuple{Int64, Combination{N}, Vector{Tuple{Int64, BernsteinTerm, Int, Combination{N}}}}[]
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      dF_bubbles = _PΛ_F_bubble_indices(r,k,D,F,i)
      push!(d_F_bubbles, (d, F, dF_bubbles))
      i += length(dF_bubbles)
    end
  end
  @check i == binomial(r+k,k)*binomial(D+r,D-k)
  :( $(d_F_bubbles) )
end

P⁻Λ_bubble_indices(r,k,D) = P⁻Λ_bubble_indices(Val(r),Val(k),Val(D))
@generated function P⁻Λ_bubble_indices(::Val{r},::Val{k},::Val{D}) where {r,k,D}
  i=0
  N = D + 1
  d_F_bubbles = Tuple{Int64, Combination{N}, Vector{Tuple{Int64, BernsteinTerm, Int, Combination{N}}}}[]
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      dF_bubbles = _P⁻Λ_F_bubble_indices(r,k,D,F,i)
      push!(d_F_bubbles, (d, F, dF_bubbles))
      i += length(dF_bubbles)
    end
  end
  @check i == binomial(r+k-1,k)*binomial(D+r,D-k)
  :( $(d_F_bubbles) )
end

# TODO Apply Hodge
sorted_and_sub_combinations(::Val,::Val{0}) = Tuple{}[]
@generated function sorted_and_sub_combinations(::Val{D},::Val{k}) where {D,k}
  @check k>0
  res = Tuple{Int, Combination{D}, Vector{Tuple{Int, Combination{D}, Int}}}[]
  for (I_id,I) in enumerate(sorted_combinations(D,k))
    sub_combis = Tuple{Int, Combination{D}, Int}[]
    for (q,I_sub_q) in enumerate(sub_combinations(I))
      I_sub_q_id = combination_index(I_sub_q)
      push!(sub_combis, (q, I_sub_q, I_sub_q_id) )
    end
    push!(res, (I_id, I, sub_combis))
  end
  :( $(res) )
end

@generated function hodged_basis_forms(::Val{D},::Val{k}) where {D,k}
  res = Tuple{Int,Int,Combination{D}}[]
  for (_,I) in enumerate(sorted_combinations(Val(D),Val(k)))
    Icomp = complement(I)
    Istar_id = combination_index(Icomp)
    sgnIcomp = combination_sign(I)
    # actually ⋆(I) = Istar = sgnIcomp*Icomp
    push!(res, (sgnIcomp, Istar_id, I) )
  end
  :( $(res) )
end
