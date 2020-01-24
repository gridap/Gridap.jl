
"""
    struct QGradMonomialBasis{...} <: Field

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Nedelec reference elements on n-cubes.
The type parameters and fields of this `struct` are not public.  
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct QGradMonomialBasis{D,T} <: Field
  order::Int
  terms::CartesianIndices{D}
  perms::Matrix{Int}
  function QGradMonomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    new{D,T}(order,terms,perms)
  end
end

"""
    QGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}

Returns a `QGradMonomialBasis` object. `D` is the dimension
of the coordinate space and `T` is the type of the components in the vector-value.
"""
function QGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}
  @assert T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  _t = tfill(order+1,Val{D-1}())
  t = (order,_t...)
  terms = CartesianIndices(t)
  perms = _prepare_perms(D)
  QGradMonomialBasis(T,order,terms,perms)
end

function field_cache(f::QGradMonomialBasis{D,T},x) where {D,T}
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_qgrad(f)
  n = 1 + f.order
  V = VectorValue{D,T}
  r = CachedArray(zeros(V,(np,ndof)))
  v = CachedArray(zeros(V,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  (r, v, c)
end

function evaluate_field!(cache,f::QGradMonomialBasis{D,T},x) where {D,T}
  r, v, c = cache
  np = length(x)
  ndof = _ndofs_qgrad(f)
  n = 1 + f.order
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
      _evaluate_nd_qgrad!(v,xi,f.order,f.terms,f.perms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r
end

function gradient_cache(f::QGradMonomialBasis{D,T},x) where {D,T}
  @assert D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_qgrad(f)
  n = 1 + f.order
  xi = testitem(x)
  V = VectorValue{D,T}
  G = gradient_type(V,xi)
  r = CachedArray(zeros(G,(np,ndof)))
  v = CachedArray(zeros(G,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  g = CachedArray(zeros(T,(D,n)))
  (r, v, c, g)
end

function evaluate_gradient!(cache,f::QGradMonomialBasis{D,T},x) where {D,T}
  r, v, c, g = cache
  np = length(x)
  ndof = _ndofs_qgrad(f)
  n = 1 + f.order
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  V = VectorValue{D,T}
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_qgrad!(v,xi,f.order,f.terms,f.perms,c,g,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r
end

# Helpers

#_ndofs_qgrad(f::QGradMonomialBasis{D}) where D = D*f.order*(f.order+1)^(D-1)

_ndofs_qgrad(f::QGradMonomialBasis{D}) where D = D*(length(f.terms))

function _prepare_perms(D)
  perms = zeros(Int,D,D)
  for j in 1:D
    for d in j:D
      perms[d,j] =  d-j+1
    end
    for d in 1:(j-1)
      perms[d,j] =  d+(D-j)+1
    end
  end
  perms
end

function _evaluate_nd_qgrad!(
  v::AbstractVector{V},
  x,
  order,
  terms::CartesianIndices{D},
  perms::Matrix{Int},
  c::AbstractMatrix{T}) where {V,T,D}

  dim = D
  for d in 1:dim
    _evaluate_1d!(c,x,order,d)
  end

  o = one(T)
  k = 1
  m = zero(mutable(V))
  js = eachindex(m)
  z = zero(T)

  for ci in terms

    for j in js

      @inbounds for i in js
        m[i] = z
      end

      s = o
      @inbounds for d in 1:dim
        s *= c[d,ci[perms[d,j]]]
      end

      m[j] = s
      v[k] = m
      k += 1

    end

  end

end

function _gradient_nd_qgrad!(
  v::AbstractVector{G},
  x,
  order,
  terms::CartesianIndices{D},
  perms::Matrix{Int},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    _evaluate_1d!(c,x,order,d)
    _gradient_1d!(g,x,order,d)
  end

  z = zero(mutable(V))
  m = zero(mutable(G))
  js = eachindex(z)
  mjs = eachindex(m)
  o = one(T)
  zi = zero(T)
  k = 1

  for ci in terms

    for j in js

      s = z
      for i in js
        s[i] = o
      end

      for q in 1:dim
        for d in 1:dim
          if d != q
            @inbounds s[q] *= c[d,ci[perms[d,j]]]
          else
            @inbounds s[q] *= g[d,ci[perms[d,j]]]
          end
        end
      end

      @inbounds for i in mjs
        m[i] = zi
      end

      for i in js
        @inbounds m[i,j] = s[i]
      end
        @inbounds v[k] = m
      k += 1

    end

  end

end
