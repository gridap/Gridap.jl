struct NedelecPrebasisOnSimplex{D} <: AbstractVector{Monomial}
  order::Int
  function NedelecPrebasisOnSimplex{D}(order::Integer) where D
    new{D}(Int(order))
  end
end

function Base.size(a::NedelecPrebasisOnSimplex{D}) where D
  k = a.order+1
  n = div(k*prod(i->(k+i),2:D),factorial(D-1))
  (n,)
end

Base.getindex(a::NedelecPrebasisOnSimplex,i::Integer) = Monomial()
Base.IndexStyle(::Type{<:NedelecPrebasisOnSimplex}) = IndexLinear()

num_terms(a::NedelecPrebasisOnSimplex) = length(a)
get_order(f::NedelecPrebasisOnSimplex) = f.order

return_type(::NedelecPrebasisOnSimplex{D}) where {D} = VectorValue{D,Float64}

function return_cache(
  f::NedelecPrebasisOnSimplex{D},x::AbstractVector{<:Point}) where D
  np = length(x)
  ndofs = num_terms(f)
  V = eltype(x)
  a = zeros(V,(np,ndofs))
  k = f.order+1
  P = MonomialBasis(Val(D),VectorValue{D,Float64},k-1,(e,order)->sum(e)<=order)
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
  mb = MonomialBasis(Val(D),VectorValue{D,Float64},k-1,(e,order)->sum(e)<=order)
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

