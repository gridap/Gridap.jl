"""
    NedelecPolyBasisOnSimplex{D,V,K,PT} <: PolynomialBasis{D,V,K,PT}

Basis of the vector valued (`V<:VectorValue{D}`) space ℕ𝔻ᴰₙ(△) for `D`=2,3.
This space is the polynomial space for Nedelec elements on simplices with
curl in (ℙᴰₙ)ᴰ. Its maximum degree is n+1 = `K`. `get_order` on it returns `K`.

Currently, the basis is implemented as the union of a UniformPolyBasis{...,PT}
for ℙᴰₙ and a monomial basis for x × (ℙᴰₙ \\ ℙᴰₙ₋₁)ᴰ.
"""
struct NedelecPolyBasisOnSimplex{D,V,K,PT} <: PolynomialBasis{D,V,K,PT}
  order::Int
  function NedelecPolyBasisOnSimplex{D}(::Type{PT},::Type{T},order::Integer) where {D,PT<:Polynomial,T}
    @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
    @check isconcretetype(PT) "PT needs to be a concrete <:Polynomial type"
    @notimplementedif !(D in (2,3))
    K = Int(order)+1
    V = VectorValue{D,T}
    new{D,V,K,PT}(Int(order))
  end
end

function Base.size(a::NedelecPolyBasisOnSimplex{D}) where D
  k = a.order+1
  n = div(k*prod(i->(k+i),2:D),factorial(D-1))
  (n,)
end

function return_cache(
  f::NedelecPolyBasisOnSimplex{D,V,K,PT},x::AbstractVector{<:Point}) where {D,V,K,PT}

  np = length(x)
  ndofs = length(f)
  a = zeros(V,(np,ndofs))
  P = UniformPolyBasis(PT,Val(D),V,K-1,_p_filter)
  cP = return_cache(P,x)
  CachedArray(a), cP, P
end

function evaluate!(
  cache,f::NedelecPolyBasisOnSimplex{3,V},x::AbstractVector{<:Point}) where V
  ca,cP,P = cache
  K = get_order(f)
  np = length(x)
  ndofs = length(f)
  ndofsP = length(P)
  setsize!(ca,(np,ndofs))
  Px = evaluate!(cP,P,x)
  a = ca.array
  T = eltype(V)
  z = zero(T)
  for (i,xi) in enumerate(x)
    # terms for (ℙₖ)³
    for j in 1:ndofsP
      a[i,j] = Px[i,j]
    end
    # terms for x × (ℙₖ\ℙₖ₋₁)³
    j = ndofsP
    x1,x2,x3 = x[i]
    for β in 1:K
      for α in 1:(K+1-β)
        j += 1
        a[i,j] = VectorValue(
          -x1^(α-1)*x2^(K-α-β+2)*x3^(β-1),
           x1^α*x2^(K-α-β+1)*x3^(β-1),
           z)
        j += 1
        a[i,j] = VectorValue(
          -x1^(K-α-β+1)*x2^(β-1)*x3^α,
          z,
          x1^(K-α-β+2)*x2^(β-1)*x3^(α-1))
      end
    end
    for γ in 1:K
      j += 1
      a[i,j] = VectorValue(
        z,
        -x2^(γ-1)*x3^(K-γ+1),
        x2^γ*x3^(K-γ))
    end
  end
  a
end

function evaluate!(
  cache,f::NedelecPolyBasisOnSimplex{2,V},x::AbstractVector{<:Point}) where V
  ca,cP,P = cache
  K = get_order(f)
  np = length(x)
  ndofs = length(f)
  ndofsP = length(P)
  setsize!(ca,(np,ndofs))
  a = ca.array
  T = eltype(V)
  z = zero(T)
  Px = evaluate!(cP,P,x)
  for (i,xi) in enumerate(x)
    # terms for (ℙₖ)²
    for j in 1:ndofsP
      a[i,j] = Px[i,j]
    end
    # terms for x × (ℙₖ\ℙₖ₋₁)²
    j = ndofsP
    x1,x2 = xi
    for α in 1:K
      j += 1
      a[i,j] = VectorValue(-x1^(α-1)*x2^(K-α+1),x1^α*x2^(K-α))
    end
    #u = one(T)
    #a[i,1] = VectorValue((u,z))
    #a[i,2] = VectorValue((z,u))
    #a[i,3] = VectorValue((-xi[2],xi[1]))
  end
  a
end

function return_cache(
  g::FieldGradientArray{1,<:NedelecPolyBasisOnSimplex{D,V,K,PT}},
  x::AbstractVector{<:Point}) where {D,V,K,PT}
  f = g.fa
  np = length(x)
  ndofs = length(f)
  xi = testitem(x)
  G = gradient_type(V,xi)
  a = zeros(G,(np,ndofs))
  mb = UniformPolyBasis(PT,Val(D),V,K-1,_p_filter)
  P = Broadcasting(∇)(mb)
  cP = return_cache(P,x)
  CachedArray(a), cP, P
end

function evaluate!(
  cache,
  g::FieldGradientArray{1,<:NedelecPolyBasisOnSimplex{3,V}},
  x::AbstractVector{<:Point}) where V
  ca,cP,P = cache
  f = g.fa
  K = get_order(f)
  np = length(x)
  ndofs = length(f)
  setsize!(ca,(np,ndofs))
  a = ca.array
  fill!(a,zero(eltype(a)))
  ndofsP = length(P)
  Px = evaluate!(cP,P,x)
  T = eltype(V)
  z = zero(T)
  for (i,xi) in enumerate(x)
    # terms for ∇((ℙₖ)³)
    for j in 1:ndofsP
      a[i,j] = Px[i,j]
    end
    # terms for  ∇(x × (ℙₖ\ℙₖ₋₁)³)
    j = ndofsP
    x1,x2,x3 = x[i]
    for β in 1:K
      for α in 1:(K+1-β)
        j += 1
        a[i,j] = TensorValue(
          #-x1^(α-1)*x2^(K-α-β+2)*x3^(β-1),
          -(α-1)*_exp(x1,α-2)*x2^(K-α-β+2)*x3^(β-1),
          -x1^(α-1)*(K-α-β+2)*_exp(x2,K-α-β+1)*x3^(β-1),
          -x1^(α-1)*x2^(K-α-β+2)*(β-1)*_exp(x3,β-2),
           #x1^α*x2^(K-α-β+1)*x3^(β-1),
           α*_exp(x1,α-1)*x2^(K-α-β+1)*x3^(β-1),
           x1^α*(K-α-β+1)*_exp(x2,K-α-β)*x3^(β-1),
           x1^α*x2^(K-α-β+1)*(β-1)*_exp(x3,β-2),
           z,z,z)
        j += 1
        a[i,j] = TensorValue(
          #-x1^(K-α-β+1)*x2^(β-1)*x3^α,
          -(K-α-β+1)*_exp(x1,K-α-β)*x2^(β-1)*x3^α,
          -x1^(K-α-β+1)*(β-1)*_exp(x2,β-2)*x3^α,
          -x1^(K-α-β+1)*x2^(β-1)*α*_exp(x3,α-1),
          z,z,z,
          #x1^(K-α-β+2)*x2^(β-1)*x3^(α-1),
          (K-α-β+2)*_exp(x1,K-α-β+1)*x2^(β-1)*x3^(α-1),
          x1^(K-α-β+2)*(β-1)*_exp(x2,β-2)*x3^(α-1),
          x1^(K-α-β+2)*x2^(β-1)*(α-1)*_exp(x3,α-2))
      end
    end
    for γ in 1:K
      j += 1
      a[i,j] = TensorValue(
        z,z,z,
        #-x2^(γ-1)*x3^(K-γ+1),
        -0*x2^(γ-1)*x3^(K-γ+1),
        -(γ-1)*_exp(x2,γ-2)*x3^(K-γ+1),
        -x2^(γ-1)*(K-γ+1)*_exp(x3,K-γ),
        #x2^γ*x3^(K-γ),
        0*x2^γ*x3^(K-γ),
        γ*_exp(x2,γ-1)*x3^(K-γ),
        x2^γ*(K-γ)*_exp(x3,K-γ-1))
    end
    #u = one(T)
    #a[i,4] = TensorValue((z,-u,z, u,z,z, z,z,z))
    #a[i,5] = TensorValue((z,z,-u, z,z,z, u,z,z))
    #a[i,6] = TensorValue((z,z,z, z,z,-u, z,u,z))
  end
  a
end

_exp(a,y) = y>0 ? a^y : one(a)

function evaluate!(
  cache,
  g::FieldGradientArray{1,<:NedelecPolyBasisOnSimplex{2,V}},
  x::AbstractVector{<:Point}) where V
  f = g.fa
  ca,cP,P = cache
  K = get_order(f)
  np = length(x)
  ndofs = length(f)
  setsize!(ca,(np,ndofs))
  a = ca.array
  fill!(a,zero(eltype(a)))
  T = eltype(V)
  z = zero(T)
  ndofsP = length(P)
  Px = evaluate!(cP,P,x)
  for (i,xi) in enumerate(x)
    # terms for  ∇((ℙₖ)²)
    for j in 1:ndofsP
      a[i,j] = Px[i,j]
    end
    # terms for  ∇(x × (ℙₖ\ℙₖ₋₁)²)
    j = ndofsP
    x1,x2 = x[i]
    for α in 1:K
      j += 1
      a[i,j] = TensorValue(
        #-x1^(α-1)*x2^(K-α+1),
        -(α-1)*_exp(x1,α-2)*x2^(K-α+1),
        -x1^(α-1)*(K-α+1)*_exp(x2,K-α),
        #x1^α*x2^(K-α),
        α*_exp(x1,α-1)*x2^(K-α),
        x1^α*(K-α)*_exp(x2,K-α-1))
    end
    #u = one(T)
    #a[i,3] = TensorValue((z,-u, u,z))
  end
  a
end

####################################
# Basis for Nedelec on D-simplices #
####################################

"""
    PGradBasis(::Type{Monomial}, ::Val{D}, ::Type{T}, order::Int) :: PolynomialBasis

Return a basis of

ℕ𝔻ᴰₙ(△) = (ℙᴰₙ)ᴰ ⊕ x × (ℙᴰₙ \\ ℙᴰₙ₋₁)ᴰ

with n=`order`, the polynomial space for Nedelec elements on `D`-dimensional
simplices with scalar type `T`. `D` must be 1, 2 or 3.

The `order`=n argument has the following meaning: the curl of the  functions in
this basis is in (ℙᴰₙ)ᴰ.

# Example:

```jldoctest
# a basis for Nedelec on tetrahedra with curl in ℙ₂
b = PGradBasis(Monomial, Val(3), Float64, 2)
```
"""
function PGradBasis(::Type{PT},::Val{D},::Type{T},order::Int) where {PT,D,T}
  # Although NedelecPolyBasisOnSimplex can be constructed with any PT<Polynomial,
  # the code explicitely uses monomials for the terms of  x×(ℙₙ \\ ℙₙ₋₁)ᴰ, so I
  # disable them here.
  # But one can use NedelecPolyBasisOnSimplex{D}(PT,T,order) if they wish.
  @notimplemented "Nedelec on simplices is only implemented for monomials"
end

