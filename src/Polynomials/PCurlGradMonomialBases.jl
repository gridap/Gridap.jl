
"""
struct PCurlGradMonomialBasis{...} <: AbstractArray{Monomial}

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Raviart-Thomas reference elements on simplices.
The type parameters and fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct PCurlGradMonomialBasis{D,T} <: AbstractVector{Monomial}
  order::Int
  pterms::Array{CartesianIndex{D},1}
  sterms::Array{CartesianIndex{D},1}
  perms::Matrix{Int}
  function PCurlGradMonomialBasis(::Type{T},order::Int,
      pterms::Array{CartesianIndex{D},1},sterms::Array{CartesianIndex{D},1},
      perms::Matrix{Int}) where {D,T}
    new{D,T}(order,pterms,sterms,perms)
  end
end

Base.size(a::PCurlGradMonomialBasis) = (_ndofs_pgrad(a),)
# @santiagobadia : Not sure we want to create the monomial machinery
Base.getindex(a::PCurlGradMonomialBasis,i::Integer) = Monomial()
Base.IndexStyle(::PCurlGradMonomialBasis) = IndexLinear()

"""
PCurlGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}

Returns a `PCurlGradMonomialBasis` object. `D` is the dimension
of the coordinate space and `T` is the type of the components in the vector-value.
The `order` argument has the following meaning: the divergence of the  functions
in this basis is in the P space of degree `order`.
"""
function PCurlGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  P_k = MonomialBasis{D}(T, order, _p_filter)
  S_k = MonomialBasis{D}(T, order, _s_filter)
  pterms = P_k.terms
  sterms = S_k.terms
  perms = _prepare_perms(D)
  PCurlGradMonomialBasis(T,order,pterms,sterms,perms)
end

"""
    num_terms(f::PCurlGradMonomialBasis{D,T}) where {D,T}
"""
function num_terms(f::PCurlGradMonomialBasis{D,T}) where {D,T}
  Int(_p_dim(f.order,D)*D + _p_dim(f.order,D-1))
end

get_order(f::PCurlGradMonomialBasis{D,T}) where {D,T} = f.order

return_type(::PCurlGradMonomialBasis{D,T}) where {D,T} = T

function return_cache(f::PCurlGradMonomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_pgrad(f)
  n = 1 + f.order+1
  V = VectorValue{D,T}
  r = CachedArray(zeros(V,(np,ndof)))
  v = CachedArray(zeros(V,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  (r, v, c)
end

function evaluate!(cache,f::PCurlGradMonomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  r, v, c = cache
  np = length(x)
  ndof = _ndofs_pgrad(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
      _evaluate_nd_pcurlgrad!(v,xi,f.order+1,f.pterms,f.sterms,f.perms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function return_cache(
  fg::FieldGradientArray{1,PCurlGradMonomialBasis{D,T}},
  x::AbstractVector{<:Point})  where {D,T}

  f = fg.fa
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_pgrad(f)
  n = 1 + f.order+1
  xi = testitem(x)
  V = VectorValue{D,T}
  G = gradient_type(V,xi)
  r = CachedArray(zeros(G,(np,ndof)))
  v = CachedArray(zeros(G,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  g = CachedArray(zeros(T,(D,n)))
  (r, v, c, g)
end

function evaluate!(cache,
  fg::FieldGradientArray{1,PCurlGradMonomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  r, v, c, g = cache
  np = length(x)
  ndof = _ndofs_pgrad(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  V = VectorValue{D,T}
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_pcurlgrad!(v,xi,f.order+1,f.pterms,f.sterms,f.perms,c,g,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end


# Helpers

_p_filter(e,order) = (sum(e) <= order)
_s_filter(e,order) = (sum(e) == order)

function _p_dim(order,D)
  dim = 1
  for d in 1:D
    dim *= order+d
  end
  dim/factorial(D)
end

_ndofs_pgrad(f::PCurlGradMonomialBasis{D}) where D = num_terms(f)


function _evaluate_nd_pcurlgrad!(
  v::AbstractVector{V},
  x,
  order,
  pterms::Array{CartesianIndex{D},1},
  sterms::Array{CartesianIndex{D},1},
  perms::Matrix{Int},
  c::AbstractMatrix{T}) where {V,T,D}

  dim = D
  for d in 1:dim
    _evaluate_1d!(c,x,order,d)
  end

  o = one(T)
  k = 1
  m = zero(Mutable(V))
  js = eachindex(m)
  z = zero(T)

  for ci in pterms
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

  for ci in sterms
    @inbounds for i in js
      m[i] = z
    end
    for j in js

      s = c[j,2]
      @inbounds for d in 1:dim
        s *= c[d,ci[d]]
      end

      m[j] = s

    end
    v[k] = m
    k += 1
  end
end

function _gradient_nd_pcurlgrad!(
  v::AbstractVector{G},
  x,
  order,
  pterms::Array{CartesianIndex{D},1},
  sterms::Array{CartesianIndex{D},1},
  perms::Matrix{Int},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    _evaluate_1d!(c,x,order,d)
    _gradient_1d!(g,x,order,d)
  end

  z = zero(Mutable(V))
  m = zero(Mutable(G))
  js = eachindex(z)
  mjs = eachindex(m)
  o = one(T)
  zi = zero(T)
  k = 1

  for ci in pterms
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

  for ci in sterms

    @inbounds for i in mjs
      m[i] = zi
    end

    for j in js

      s = z
      for i in js
        s[i] = c[j,2]
      end

      for q in 1:dim
        for d in 1:dim
          if d != q
            @inbounds s[q] *= c[d,ci[d]]
          else
            @inbounds s[q] *= g[d,ci[d]]
          end
        end
      end
      aux = o
      @inbounds for d in 1:dim
        aux *= c[d,ci[d]]
      end
      s[j] += aux

      for i in js
        @inbounds m[i,j] = s[i]
      end
    end
    @inbounds v[k] = m
    k += 1
  end
end
