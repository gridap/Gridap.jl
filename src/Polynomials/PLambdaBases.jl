#################################
# Tensorial nD polynomial bases #
#################################

"""
struct PLambda{D,T,L,B,P} <: PolynomialBasis{D,VectorValue{L,T},Bernstein}

Finite Element Exterior Calculus polynomial basis for the spaces `D`-dimensional
simplices, P`r`Λ`ᴷ`(△`ᴰ`).

Reference: D. N. Arnold and A. Logg, Periodic Table of the Finite Elements, SIAM News, vol. 47 no. 9, November 2014
"""
struct PLambdaBasis{D,T,L,P,B} <: PolynomialBasis{D,VectorValue{L,T},Bernstein}
  r::Int
  k::Int
  Ψ::P
  scalar_bernstein_basis::B

  function PLambdaBasis{D}(::Type{T},r,k,vertices=nothing) where {D,T}
    @check T<:Real "T needs to be <:Real since represents the scalar type"
    @check k in 0:D "The form order k must be in 0:D"
    @check r > 0    "The polynomial order r must be positive"
    if !isnothing(vertices)
      @check length(vertices) == D+1 "$D+1 vertices are required to define a $D-dim simplex, got $(length(vertices))"
      @check eltype(vertices) isa Point{D} "Vertices should be of type Point{$D}, got $(eltype(vertices))"
    end

    L = binomial(k,D) # Number of components of a basis form
    C = binomial(r+k,k)*binomial(D+r,D-k) # Number of basis polynomials

    b = BernsteinBasisOnSimplex{D}(T, r, vertices=vertices)
    B = typeof(b)

    if isnothing(vertices)
      Ψ = nothing
      P = Nothing
    else
      # Col major: the L components of a basis form are stored contiguously
      P = SMatrix{L,C,T}
      Ψ = MMatrix{L,C,T}(undef)
      Ψ = _compute_PΛ_basis_form_coefficient!(Ψ,r,k,D,b,vertices)
    end

    new{D,T,L,P,B}(r,k,Ψ,b)
  end
end

function PLambdaBasis(::Val{D},::Type{T},r,k,vertices=nothing) where {D,T}
  PLambdaBasis{D}(T,r,vertices)
end

get_FEEC_poly_degree(b::PLambdaBasis) = b.r
get_FEEC_form_degree(b::PLambdaBasis) = b.k
get_FEEC_family(::PLambdaBasis) = :P

function Base.size(b::PLambdaBasis{D}) where D
  r, k = b.r, b.k
  return (binomial(r+k,k)*binomial(D+r,D-k), )
end
get_order(b::PLambdaBasis) = get_FEEC_poly_degree(b)

get_cart_to_bary_matrix(b::PLambdaBasis) = b.scalar_bernstein_basis.cart_to_bary_matrix
#get_dimension(::PLambdaBasis{D}) where D = D


###################
# Implementation  #
###################

_compute_PΛ_basis_form_coefficient!(Ψ,r,k,D,b,::Nothing) = nothing # done in _hat_Ψ

function _compute_PΛ_basis_form_coefficient!(Ψ,r,k,D,b,vertices)
  N = D+1
  T = eltype(Ψ)
  α_prec = ntuple(_->-1, N)
  φ_αF = MMatrix{D,N,T}(undef)
  for (_, F, dF_bubbles) in PΛ_bubble_indices(r,k,D)
    for (w, α, J) in dF_bubbles
      if α ≠ α_prec
        update_φ_αF!(φ_αF,b,α,F,r)
        α_prec = α
      end
      for (iI,I) in enumerate(sorted_combinations(Val{k},Val{D}))
        Ψ[iI,w] = @inline minor(φ_αF,I,J)
      end
    end
  end
end

@inline function update_φ_αF!(φ_αF,b,α,F,r)
  M = b.cart_to_bary_matrix
  @inbounds for ci in CartesianIndices(φ_αF)
    i, j = ci[1], ci[2]
    mF = sum(M[Fl,i+1] for Fl in F; init=0)
    φ_αF[i][j] = M[j,i+1] - α[j]*mF/r
  end
end

"""
    _hat_Ψ(r,α,F,I,J,T)::T

PLambdaBasis.Ψ matrix elements in the reference simplex, T is the scalar return type
"""
function _hat_Ψ(r,α,F,I,J,::Type{T})::T where T
  @check sum(α) == r
  @check length(I) == length(J) # thanks to dispatch on zero length J below
  @check length(J) > 0 # thanks to dispatch on zero length J below

  @inbounds begin

    s = Int(isone(J[1]))
    n = count(i-> (J[i]-1)∉I, (s+1):length(J))

    n > 1 && return 0. # rank M_IJ inferior to 2

    p = _find_first_val_or_zero(j-> (I[j]+1)∉J, 1, length(J))

    if isone(n)        # rank M_IJ is 1
      m = _find_first_val_or_zero(i-> (J[i]-1)∉I, (s+1), length(J))
      u_m, v_p = _u(m,F,I), _v(p,α,J,r)
      iszero(s) && return (-1)^(m+p)*u_m*v_p

      q = _find_first_val_or_zero(j-> (I[j]+s)∉J, (p+1), length(J))
      @check !iszero(q)
      v_q = _v(q,α,J,r)
      return (-1)^(m+p+q) * u_m * (v_q - v_p)
    end

    u, v = _u(F,I), _v(α,J,r)
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
_hat_Ψ(r,α,F,I,J::Combination{0},::Type{T}) where T = one(T) # 0 forms

@propagate_inbounds function _find_first_val_or_zero(pred, start, stop)
  r = findfirst(pred,start:stop)
  return isnothing(r) ? 0 : r+start-1
end

@propagate_inbounds _u(i,F,I)   = Int(isone(F[1])) - Int(I[i]+1 in F)
@propagate_inbounds _u(F,I::Combination{k}) where k = ntuple(i->_u(i,F,I), Val(k))
@propagate_inbounds _v(j,α,J,r) = α[J[j]]/r
@propagate_inbounds _v(α,J::Combination{k},r) where k = ntuple(j->_v(j,α,J,r), Val(k))


# API

function _return_cache(
  b::PLambdaBasis{D},x,::Type{G},::Val{N_deriv}) where {D,G,N_deriv}

  @assert D == length(eltype(x)) "Incorrect number of point components"
  T = eltype(G)
  r = get_order(b)
  np = length(x)
  ndof = length(b)
  ndof_bernstein = binomial(r+D,D)
  ndof_val = binomial(r+D,D)

  r = CachedArray(zeros(G,(np,ndof)))
  s = MArray{Tuple{Vararg{D,N_deriv}},T}(undef)
  # Cache for all scalar nD-Bernstein polynomials, no other caches needed for derivatives
  c = CachedVector(zeros(T,ndof_scalar))
  # Cache for 1
  t = ntuple( _ -> nothing, Val(N_deriv))
  (r, s, c, t...)
end

function return_cache(
  fg::FieldGradientArray{1,<:PLambdaBasis},
  x::AbstractVector{<:Point})

  fg_b = FieldGradientArray{1}(fg.fa._basis)
  return (fg_b, return_cache(fg_b,x))
end

function return_cache(
  fg::FieldGradientArray{2,<:PLambdaBasis},
  x::AbstractVector{<:Point})

  fg_b = FieldGradientArray{2}(fg.fa._basis)
  return (fg_b, return_cache(fg_b,x))
end

function _evaluate_nd!(
  b::PLambdaBasis{D,T}, x,
  ω::AbstractMatrix{V}, i,
  c::AbstractVector{T}) where {D,V,T}

  r = get_FEEC_poly_degree(b)
  k = get_FEEC_form_degree(b)
  λ = _cart_to_bary(x, get_cart_to_bary_matrix(b))

  # _evaluate_nd!(::BernsteinBasisOnSimplex) without set_value
  c[1] = one(T)
  _downwards_de_Casteljau_nD!(c,λ,Val(r),Val(D))

  for (_, _, dF_bubbles) in PΛ_bubble_indices(r,k,D)
    for (w, α, _) in dF_bubbles
      id_α = _simplex_multi_id_to_linear_id(α)
      ω[i][w] = c[id_α] .* b.Ψ[w] # Bα(x)*Ψ_w
    end
  end
end

function _gradient_nd!(
  b::PLambdaBasis{D,T}, x,
  ∇ω::AbstractMatrix{G}, i,
  c::AbstractVector{T},
  ∇B::MMatrix{1,LB,VectorValue{D,T}},
  s::MVector{D,T}) where {D,G,T,LB}

  #r = get_FEEC_poly_degree(b)
  #k = get_FEEC_form_degree(b)
  #x_to_λ = get_cart_to_bary_matrix(b)
  #λ = _cart_to_bary(x, x_to_λ)
  #c[1] = one(T)
  #_downwards_de_Casteljau_nD!(c,λ,Val(r-1),Val(D))
  #_grad_Bα_from_Bα⁻!(∇B,1,c,s,Val(r),Val(D),T,x_to_λ)

  # _gradient_nd!(::BernsteinBasisOnSimplex) with set_value in our cache ∇B
  _gradient_nd!(b.B, x, ∇B, 1, c, nothing, s)

  for (_, _, dF_bubbles) in PΛ_bubble_indices(r,k,D)
    for (w, α, _) in dF_bubbles
      id_α = _simplex_multi_id_to_linear_id(α) # TODO optimize this
      ∇ω[i][w] = ∇B[1,id_α] ⊗  VectorValue(b.Ψ[w])
    end
  end
end

function _hessian_nd!(
  b::PLambdaBasis{D,T}, x,
  Hω::AbstractMatrix{G}, i,
  c::AbstractVector{T},
  HB::MMatrix{1,LB,TensorValue{D,D,T}},
  s::MMatrix{D,T}) where {D,G,T,LB}

  #r = get_FEEC_poly_degree(b)
  #k = get_FEEC_form_degree(b)
  #x_to_λ = get_cart_to_bary_matrix(b)
  #λ = _cart_to_bary(x, x_to_λ)
  #c[1] = one(T)
  #_downwards_de_Casteljau_nD!(c,λ,Val(r-2),Val(D))
  #_hess_Bα_from_Bα⁻⁻!(HB,1,c,s,Val(r),Val(D),T,x_to_λ)

  # _hessian_nd!(::BernsteinBasisOnSimplex) with set_value in our cache HB
  _hessian_nd!(b.B, x, HB, 1, c, nothing, nothing, s)

  for (_, _, dF_bubbles) in PΛ_bubble_indices(r,k,D)
    for (w, α, _) in dF_bubbles
      id_α = _simplex_multi_id_to_linear_id(α) # TODO optimize this
      Hω[i][w] = HB[1,id_α] ⊗  VectorValue(b.Ψ[w])
    end
  end
end


##########################
# PLambda bases helpers  #
##########################


PΛ_bubble_indices(r,k,D) = PΛ_bubble_indices(Val(r),Val(k),Val(D))
@generated function PΛ_bubble_indices(::Val{r},::Val{k},::Val{D}) where {r,k,D}
  i=0
  d_F_bubbles = []
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      dF_bubbles = _PΛ_F_bubble_indices(r,k,D,F,i)
      push!(d_F_bubbles, (d, F, dF_bubbles))
      i += length(dF_bubbles)
    end
  end
  @check i == binomial(r+k,k)*binomial(D+r,D-k)
  d_F_bubbles = tuple(d_F_bubbles...)
  :( $(d_F_bubbles) )
end

P⁻Λ_bubble_indices(r,k,D) = P⁻Λ_bubble_indices(Val(r),Val(k),Val(D))
@generated function P⁻Λ_bubble_indices(::Val{r},::Val{k},::Val{D}) where {r,k,D}
  i=0
  d_F_bubbles = []
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      dF_bubbles = _P⁻Λ_F_bubble_indices(r,k,D,F,i)
      push!(d_F_bubbles, (d, F, dF_bubbles))
      i += length(dF_bubbles)
    end
  end
  @check i == binomial(r+k-1,k)*binomial(D+r,D-k)
  d_F_bubbles = tuple(d_F_bubbles...)
  :( $(d_F_bubbles) )
end

function _P⁻Λ_F_bubble_indices(r,k,D,F,i)
  N = D + 1
  ids = []
  for α in bernstein_terms(r-1,N)
    for J in sorted_combinations(N,k+1)
      j = minimum(J)-1
      if issetequal(supp(α) ∪ J, F) && all(α[1:j] .== 0)
        push!(ids, (i, α, J))
        i += 1
      end
    end
  end
  tuple(ids...)
end

function _PΛ_F_bubble_indices(r,k,D,F,i)
  N = D + 1
  ids = []
  for α in bernstein_terms(r,N)
    for J in sorted_combinations(N,k)
      j = minimum(setdiff(F,J))-1
      if issetequal(supp(α) ∪ J, F) && all(α[1:j] .== 0)
        push!(ids, (i, α, J))
        i += 1
      end
    end
  end
  tuple(ids...)
end


"""
supp(α)

TBW
"""
function supp(α)
  s = Int[]
  for (i,αi) in enumerate(α)
    if αi > 0
      push!(s, i)
    end
  end
  Tuple(s)
end

function minor(M,I,J)
  @check length(I) == length(J)
  @check I ⊆ axes(M)[1]
  @check J ⊆ axes(M)[2]

  k = length(I)
  T = eltype(M)
  m = MMatrix{k,k,T}(undef)
  for (i, Ii) in enumerate(I)
    for (j, Jj) in enumerate(J)
      @inbounds m[i,j] = M[Ii,Jj]
    end
  end
  det(m)
end

function all_k_minors!(m,M,::Val{k}) where {k}
  D = size(M)[1]
  Λᵏᴰ = sorted_combinations(Val(D),Val(k))
  @inbounds begin
    for (i, I) in enumerate(Λᵏᴰ)
      for (j, J) in enumerate(Λᵏᴰ)
        m[i,j] = @inline minor(M,I,J)
      end
    end
  end
end

