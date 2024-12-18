"""
    NonTensorRTPolyBasis{D,V,K,PT} <: PolynomialBasis{D,V,K,PT}

Basis of ℝ𝕋ᴰₙ = (𝕊ₙ)ᴰ ⊕ x (𝕊ₙ\\𝕊₍ₙ₋₁₎)
where 𝕊ₙ is a multivariate scalar polynomial space of maximum degree n = `K`-1.

ℝ𝕋ᴰₙ is the space needed for Raviart-Thomas elements with divergence in 𝕊ₙ.
Its maximum degree is n+1 = `K`. `get_order` on it returns `K`.

The multivariate scalar space 𝕊ₙ, typically ℙₙ, doesn't need to have a tensor
product structure of 1D scalar spaces. Thus, the ℝ𝕋ᴰₙ component's scalar spaces
are not tensor products either.

𝕊ₙ must admit a basis computable using products of 1D polynomials of the `PT`
type, `PT` thus needs to be hierarchical, see [`isHierarchical`](@ref).

Indeed, 𝕊ₙ is defined like a scalar valued `TensorPolynomialBasis` via the
`_filter` argument of the constructor, by default `_p_filter`.
"""
struct NonTensorRTPolyBasis{D,V,K,PT} <: PolynomialBasis{D,V,K,PT}
  pterms::Vector{CartesianIndex{D}}
  sterms::Vector{CartesianIndex{D}}

  """
      NonTensorRTPolyBasis{D}(::Type{PT}, ::Type{T}, order::Int, _filter::Function=_p_filter)

  Where `_filter` defines 𝕊ₙ and `order` = n = K-1 (cf. struct docstring).
  """
  function NonTensorRTPolyBasis{D}(
    ::Type{PT}, ::Type{T}, order::Int,
    _filter::Function=_p_filter
    ) where {PT,D,T}

    @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
    @check D > 1
    @check isHierarchical(PT) "The polynomial basis must be hierarchichal for this space."

    V = VectorValue{D,T}

    # terms defining 𝕊_order
    P_k = MonomialBasis(Val(D), T, order, _filter)
    pterms = copy(P_k.terms)
    msg =  "Some term defining `𝕊ₙ` contain a higher index than the maximum,
      `order`+1, please fix the `_filter` argument"
    @check all( pterm -> (maximum(Tuple(pterm), init=0) <= order+1), pterms) msg

    # terms defining 𝕊_order\𝕊_{order-1}
    _minus_one_order_filter = term -> _filter(Tuple(term) .- 1, order-1) # .-1 for index / degree shift
    sterms = filter(!_minus_one_order_filter, pterms)

    new{D,V,order+1,PT}(pterms,sterms)
  end
end

Base.size(a::NonTensorRTPolyBasis{D}) where {D} = D*length(a.pterms) + length(a.sterms)
Base.getindex(a::NonTensorRTPolyBasis,i::Integer) = Monomial()
Base.IndexStyle(::NonTensorRTPolyBasis) = IndexLinear()

function evaluate!(
  cache,f::NonTensorRTPolyBasis{D,V,K,PT},
  x::AbstractVector{<:Point}) where {D,V,K,PT}

  r, v, c = cache
  np = length(x)
  ndof = length(f)
  ndof_1d = get_order(f) + 1 # K+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,ndof_1d))
  for i in 1:np
    @inbounds xi = x[i]
      _evaluate_nd_rt!(PT,v,xi,K,f.pterms,f.sterms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function evaluate!(cache,
  fg::FieldGradientArray{1,NonTensorRTPolyBasis{D,V,K,PT}},
  x::AbstractVector{<:Point}) where {D,V,K,PT}

  f = fg.fa
  r, v, c, g = cache
  np = length(x)
  ndof = length(f)
  ndof_1d = get_order(f) + 1 # K+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,ndof_1d))
  setsize!(g,(D,ndof_1d))
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_rt!(PT,v,xi,K,f.pterms,f.sterms,c,g,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end


# Helpers


function _evaluate_nd_rt!(
  ::Type{PT},
  v::AbstractVector{V}, # V = VectorValue{D,T}
  x,
  order,
  pterms::Vector{CartesianIndex{D}},
  sterms::Vector{CartesianIndex{D}},
  c::AbstractMatrix{T}) where {PT,V,T,D}

  for d in 1:D
    K = Val(order)
    _evaluate_1d!(PT,K,c,x,d)
  end

  m = zero(Mutable(V))
  k = 1

  @inbounds begin
    for l in 1:D
      for ci in pterms

        s = one(T)
        for d in 1:D
          s *= c[d,ci[d]]
        end

        k = _comp_wize_set_value!(v,s,k,l)
      end
    end

    for ci in sterms
      for i in 1:D
        m[i] = zero(T)
      end

      for l in 1:D

        s = x[l]
        for d in 1:D
          s *= c[d,ci[d]]
        end

        m[l] = s
      end

      v[k] = m
      k += 1
    end
  end
end

function _gradient_nd_rt!(
  ::Type{PT},
  v::AbstractVector{G},
  x,
  order,
  pterms::Vector{CartesianIndex{D}},
  sterms::Vector{CartesianIndex{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {PT,G,T,D,V}

  for d in 1:D
    K = Val(order)
    _derivatives_1d!(PT,K,(c,g),x,d)
  end

  s = zero(Mutable(V))
  m = zero(Mutable(G))
  k = 1

  @inbounds begin
    for l in 1:D
      for ci in pterms

        for i in eachindex(s)
          s[i] = one(T)
        end

        for q in 1:D
          for d in 1:D
            if d != q
               s[q] *= c[d,ci[d]]
            else
               s[q] *= g[d,ci[d]]
            end
          end
        end

        k = _comp_wize_set_gradient!(v,s,k,Val(l),V)
      end
    end

    for ci in sterms

      for i in eachindex(m)
        m[i] = zero(T)
      end

      for l in 1:D

        for i in eachindex(s)
          s[i] = x[l]
        end

        aux = one(T)
        for q in 1:D
          aux *= c[q,ci[q]]
          for d in 1:D
            if d != q
               s[q] *= c[d,ci[d]]
            else
               s[q] *= g[d,ci[d]]
            end
          end
        end
        s[l] += aux

        for i in 1:D
           m[i,l] = s[i]
        end
      end
      v[k] = m
      k += 1
    end
  end
end

"""
    PCurlGradBasis(::Type{PT}, ::Val{D}, ::Type{T}, order::Int) :: PolynomialBasis

Return a basis of ℙ_`order` ⊕ x ⋅ (ℙ_`order`\\ℙ_{`order`-1}), the polynomial
space for Raviart-Thomas elements on `D`-dimensional simplices with scalar type `T`.

The `order` argument of this function has the following meaning: the divergence
of the functions in this basis is in the ℙ space of degree `order`.

`PT<:Polynomial` is the choice of scalar 1D polynomial basis, it must be
hierarchichal, see [`isHierarchichal`](@ref).

Returns a `NonTensorRTPolyBasis{D,VectorValue{D,T},order+1,PT}` object,
or `TensorPolynomialBasis{1,VectorValue{1,T},order+1,PT}` for `D`=1.
"""
function PCurlGradBasis(::Type{PT},::Val{D},::Type{T},order::Int) where {PT,D,T}
  NonTensorRTPolyBasis{D}(PT, T, order)
end

function PCurlGradBasis(::Type{PT},::Val{1},::Type{T},order::Int) where {PT,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"

  V = VectorValue{1,T}
  TensorPolynomialBasis(PT, Val(1), V, order+1)
end
