
"""
    struct QGradMonomialBasis{...} <: AbstractVector{Monomial}

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Nedelec reference elements on n-cubes.
The type parameters and fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct QGradMonomialBasis{D,T} <: AbstractVector{Monomial}
  order::Int
  terms::CartesianIndices{D}
  perms::Matrix{Int}
  function QGradMonomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    new{D,T}(order,terms,perms)
  end
end

Base.size(a::QGradMonomialBasis) = (_ndofs_qgrad(a),)
# @santiagobadia : Not sure we want to create the monomial machinery
Base.getindex(a::QGradMonomialBasis,i::Integer) = Monomial()
Base.IndexStyle(::QGradMonomialBasis) = IndexLinear()

"""
    QGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}

Returns a `QGradMonomialBasis` object. `D` is the dimension
of the coordinate space and `T` is the type of the components in the vector-value.
The `order` argument has the following meaning: the curl of the  functions in this basis
is in the Q space of degree `order`.
"""
function QGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  _order = order + 1
  _t = tfill(_order+1,Val{D-1}())
  t = (_order,_t...)
  terms = CartesianIndices(t)
  perms = _prepare_perms(D)
  QGradMonomialBasis(T,order,terms,perms)
end

"""
    num_terms(f::QGradMonomialBasis{D,T}) where {D,T}
"""
num_terms(f::QGradMonomialBasis{D,T}) where {D,T} = length(f.terms)*D

get_order(f::QGradMonomialBasis) = f.order

function return_cache(f::QGradMonomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_qgrad(f)
  n = 1 + f.order+1
  V = VectorValue{D,T}
  r = CachedArray(zeros(V,(np,ndof)))
  v = CachedArray(zeros(V,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  (r, v, c)
end

function evaluate!(cache,f::QGradMonomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  r, v, c = cache
  np = length(x)
  ndof = _ndofs_qgrad(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
      _evaluate_nd_qgrad!(v,xi,f.order+1,f.terms,f.perms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function return_cache(
  fg::FieldGradientArray{1,QGradMonomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_qgrad(f)
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

function evaluate!(
  cache,
  fg::FieldGradientArray{1,QGradMonomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  r, v, c, g = cache
  np = length(x)
  ndof = _ndofs_qgrad(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  V = VectorValue{D,T}
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_qgrad!(v,xi,f.order+1,f.terms,f.perms,c,g,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

# Helpers

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
  m = zero(Mutable(V))
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

  z = zero(Mutable(V))
  m = zero(Mutable(G))
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

struct NedelecPrebasisOnSimplex{D} <: AbstractVector{Monomial}
  order::Int
  function NedelecPrebasisOnSimplex{D}(order::Integer) where D
    new{D}(Int(order))
  end
end

function Base.size(a::NedelecPrebasisOnSimplex{d}) where d
  k = a.order+1
  n = div(k*prod(i->(k+i),2:d),factorial(d-1))
  (n,)
end

Base.getindex(a::NedelecPrebasisOnSimplex,i::Integer) = Monomial()
Base.IndexStyle(::Type{<:NedelecPrebasisOnSimplex}) = IndexLinear()

num_terms(a::NedelecPrebasisOnSimplex) = length(a)
get_order(f::NedelecPrebasisOnSimplex) = f.order

function return_cache(
  f::NedelecPrebasisOnSimplex{d},x::AbstractVector{<:Point}) where d
  np = length(x)
  ndofs = num_terms(f)
  V = eltype(x)
  a = zeros(V,(np,ndofs))
  k = f.order+1
  P = MonomialBasis{d}(VectorValue{d,Float64},k-1,(e,order)->sum(e)<=order)
  cP = return_cache(P,x)
  CachedArray(a), cP, P
end

function evaluate!(
  cache,f::NedelecPrebasisOnSimplex{3},x::AbstractVector{<:Point})
  ca,cP,P = cache
  k = f.order+1
  np = length(x)
  ndofs = num_terms(f)
  ndofsP = length(P)
  setsize!(ca,(np,ndofs))
  Px = evaluate!(cP,P,x)
  a = ca.array
  V = eltype(x)
  T = eltype(V)
  z = zero(T)
  u = one(T)
  for (ip,p) in enumerate(x)
    for j in 1:ndofsP
      a[ip,j] = Px[ip,j]
    end
    i = ndofsP
    x1,x2,x3 = x[ip]
    zp = zero(x1)
    for β in 1:k
      for α in 1:(k+1-β)
        i += 1
        a[ip,i] = VectorValue(
          -x1^(α-1)*x2^(k-α-β+2)*x3^(β-1),
           x1^α*x2^(k-α-β+1)*x3^(β-1),
           zp)
        i += 1
        a[ip,i] = VectorValue(
          -x1^(k-α-β+1)*x2^(β-1)*x3^α,
          zp,
          x1^(k-α-β+2)*x2^(β-1)*x3^(α-1))
      end
    end
    for γ in 1:k
      i += 1
      a[ip,i] = VectorValue(
        zp,
        -x2^(γ-1)*x3^(k-γ+1),
        x2^γ*x3^(k-γ))
    end
  end
  a
end

function evaluate!(
  cache,f::NedelecPrebasisOnSimplex{2},x::AbstractVector{<:Point})
  ca,cP,P = cache
  k = f.order+1
  np = length(x)
  ndofs = num_terms(f)
  ndofsP = length(P)
  setsize!(ca,(np,ndofs))
  a = ca.array
  V = eltype(x)
  T = eltype(V)
  z = zero(T)
  u = one(T)
  Px = evaluate!(cP,P,x)
  for (ip,p) in enumerate(x)
    for j in 1:ndofsP
      a[ip,j] = Px[ip,j]
    end
    i = ndofsP
    x1,x2 = x[ip]
    zp = zero(x1)
    for α in 1:k
      i += 1
      a[ip,i] = VectorValue(-x1^(α-1)*x2^(k-α+1),x1^α*x2^(k-α))
    end
    #a[ip,1] = VectorValue((u,z))
    #a[ip,2] = VectorValue((z,u))
    #a[ip,3] = VectorValue((-p[2],p[1]))
  end
  a
end

function return_cache(
  g::FieldGradientArray{1,<:NedelecPrebasisOnSimplex{D}},
  x::AbstractVector{<:Point}) where D
  f = g.fa
  np = length(x)
  ndofs = num_terms(f)
  xi = testitem(x)
  V = eltype(x)
  G = gradient_type(V,xi)
  a = zeros(G,(np,ndofs))
  k = f.order+1
  mb = MonomialBasis{D}(VectorValue{D,Float64},k-1,(e,order)->sum(e)<=order)
  P = Broadcasting(∇)(mb)
  cP = return_cache(P,x)
  CachedArray(a), cP, P
end

function evaluate!(
  cache,
  g::FieldGradientArray{1,<:NedelecPrebasisOnSimplex{3}},
  x::AbstractVector{<:Point})
  ca,cP,P = cache
  f = g.fa
  k = f.order+1
  np = length(x)
  ndofs = num_terms(f)
  setsize!(ca,(np,ndofs))
  a = ca.array
  fill!(a,zero(eltype(a)))
  ndofsP = length(P)
  Px = evaluate!(cP,P,x)
  V = eltype(x)
  T = eltype(V)
  z = zero(T)
  u = one(T)
  for (ip,p) in enumerate(x)
    for j in 1:ndofsP
      a[ip,j] = Px[ip,j]
    end
    i = ndofsP
    x1,x2,x3 = x[ip]
    zp = zero(x1)
    for β in 1:k
      for α in 1:(k+1-β)
        i += 1
        a[ip,i] = TensorValue(
          #-x1^(α-1)*x2^(k-α-β+2)*x3^(β-1),
          -(α-1)*_exp(x1,α-2)*x2^(k-α-β+2)*x3^(β-1),
          -x1^(α-1)*(k-α-β+2)*_exp(x2,k-α-β+1)*x3^(β-1),
          -x1^(α-1)*x2^(k-α-β+2)*(β-1)*_exp(x3,β-2),
           #x1^α*x2^(k-α-β+1)*x3^(β-1),
           α*_exp(x1,α-1)*x2^(k-α-β+1)*x3^(β-1),
           x1^α*(k-α-β+1)*_exp(x2,k-α-β)*x3^(β-1),
           x1^α*x2^(k-α-β+1)*(β-1)*_exp(x3,β-2),
           #zp,
           zp,zp,zp)
        i += 1
        a[ip,i] = TensorValue(
          #-x1^(k-α-β+1)*x2^(β-1)*x3^α,
          -(k-α-β+1)*_exp(x1,k-α-β)*x2^(β-1)*x3^α,
          -x1^(k-α-β+1)*(β-1)*_exp(x2,β-2)*x3^α,
          -x1^(k-α-β+1)*x2^(β-1)*α*_exp(x3,α-1),
          # zp
          zp,zp,zp,
          #x1^(k-α-β+2)*x2^(β-1)*x3^(α-1),
          (k-α-β+2)*_exp(x1,k-α-β+1)*x2^(β-1)*x3^(α-1),
          x1^(k-α-β+2)*(β-1)*_exp(x2,β-2)*x3^(α-1),
          x1^(k-α-β+2)*x2^(β-1)*(α-1)*_exp(x3,α-2))
      end
    end
    for γ in 1:k
      i += 1
      a[ip,i] = TensorValue(
        #zp
        zp,zp,zp,
        #-x2^(γ-1)*x3^(k-γ+1),
        -0*x2^(γ-1)*x3^(k-γ+1),
        -(γ-1)*_exp(x2,γ-2)*x3^(k-γ+1),
        -x2^(γ-1)*(k-γ+1)*_exp(x3,k-γ),
        #x2^γ*x3^(k-γ),
        0*x2^γ*x3^(k-γ),
        γ*_exp(x2,γ-1)*x3^(k-γ),
        x2^γ*(k-γ)*_exp(x3,k-γ-1))
    end
    #a[ip,4] = TensorValue((z,-u,z, u,z,z, z,z,z))
    #a[ip,5] = TensorValue((z,z,-u, z,z,z, u,z,z))
    #a[ip,6] = TensorValue((z,z,z, z,z,-u, z,u,z))
  end
  a
end

_exp(a,y) = y>0 ? a^y : one(a)

function evaluate!(
  cache,
  g::FieldGradientArray{1,<:NedelecPrebasisOnSimplex{2}},
  x::AbstractVector{<:Point})
  f = g.fa
  ca,cP,P = cache
  k = f.order+1
  np = length(x)
  ndofs = num_terms(f)
  setsize!(ca,(np,ndofs))
  a = ca.array
  fill!(a,zero(eltype(a)))
  V = eltype(x)
  T = eltype(V)
  z = zero(T)
  u = one(T)
  ndofsP = length(P)
  Px = evaluate!(cP,P,x)
  for (ip,p) in enumerate(x)
    for j in 1:ndofsP
      a[ip,j] = Px[ip,j]
    end
    i = ndofsP
    x1,x2 = x[ip]
    zp = zero(x1)
    for α in 1:k
      i += 1
      a[ip,i] = TensorValue(
        #-x1^(α-1)*x2^(k-α+1),
        -(α-1)*_exp(x1,α-2)*x2^(k-α+1),
        -x1^(α-1)*(k-α+1)*_exp(x2,k-α),
        #x1^α*x2^(k-α),
        α*_exp(x1,α-1)*x2^(k-α),
        x1^α*(k-α)*_exp(x2,k-α-1))
    end
    #a[ip,3] = TensorValue((z,-u, u,z))
  end
  a
end

