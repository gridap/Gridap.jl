"""
    NedelecPolyBasisOnSimplex{D,V,PT} <: PolynomialBasis{D,V,PT}

Basis of the vector valued (`V<:VectorValue{D}`) space 𝓝𝓓ᴰₙ(△) for `D`=2,3.
This space is the polynomial space for Nedelec elements on simplices with
curl in (𝓟ᴰₙ)ᴰ. Its maximum degree is n+1 = `K`. `get_order` on it returns `K`.

   𝓝𝓓ᴰₙ(△) = (𝓟ᴰₙ)ᴰ ⊕ x × (𝓟ᴰₙ \\ 𝓟ᴰₙ₋₁)ᴰ

Currently, the basis is implemented as the union of a CartProdPolyBasis{...,PT}
for 𝓟ᴰₙ and a monomial basis for x × (𝓟ᴰₙ \\ 𝓟ᴰₙ₋₁)ᴰ.

!!! warning
    Using this basis is not recommanded, [`BarycentricPmΛBasis`](@ref) is better
    numerically conditioned for higher degrees, they are obtained by using
    `Bernstein` as argument of [`FEEC_poly_basis`](@ref) .

# Examples
These return instances of `NedelecPolyBasisOnSimplex`
```jldoctest
# a basis for Nedelec on triangles with curl in 𝓟²₁
b = FEEC_poly_basis(Val(2),Float64,2,1,:P⁻,Monomial)

# a basis for Nedelec on tetrahedra with curl in 𝓟³₁
b = FEEC_poly_basis(Val(3),Float64,2,1,:P⁻,Monomial)
```
"""
struct NedelecPolyBasisOnSimplex{D,V,PT} <: PolynomialBasis{D,V,PT}
  order::Int
  function NedelecPolyBasisOnSimplex{D}(::Type{PT},::Type{T},order::Integer) where {D,PT<:Polynomial,T}
    @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
    @check isconcretetype(PT) "PT needs to be a concrete <:Polynomial type"
    @check isHierarchical(PT) "The polynomial basis must be hierarchical for this space."
    @notimplementedif !(D in (2,3))
    V = VectorValue{D,T}
    new{D,V,PT}(Int(order))
  end
end

get_order(f::NedelecPolyBasisOnSimplex) = f.order + 1 # Return actual maximum poly order
get_orders(b::NedelecPolyBasisOnSimplex{D}) where D = tfill(get_order(b), Val(D))

function Base.size(f::NedelecPolyBasisOnSimplex{D}) where D
  K = get_order(f)
  n = div(K*prod(i->(K+i),2:D),factorial(D-1))
  (n,)
end

function testvalue(::Type{NedelecPolyBasisOnSimplex{D,V,PT}}) where {D,V,PT}
  T = eltype(V)
  NedelecPolyBasisOnSimplex{D}(PT,T,0)
end

function return_cache(
  f::NedelecPolyBasisOnSimplex{D,V,PT},x::AbstractVector{<:Point}) where {D,V,PT}

  K = get_order(f)
  np = length(x)
  ndofs = length(f)
  Vr = _return_val_eltype(f,x)
  a = zeros(Vr,(np,ndofs))
  P = CartProdPolyBasis(PT,Val(D),V,K-1,_p_filter)
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
    # terms for (𝓟ₖ)³
    for j in 1:ndofsP
      a[i,j] = Px[i,j]
    end
    # terms for x × (𝓟ₖ\𝓟ₖ₋₁)³
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
    # terms for (𝓟ₖ)²
    for j in 1:ndofsP
      a[i,j] = Px[i,j]
    end
    # terms for x × (𝓟ₖ\𝓟ₖ₋₁)²
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
  g::FieldGradientArray{1,<:NedelecPolyBasisOnSimplex{D,V,PT}},
  x::AbstractVector{<:Point}) where {D,V,PT}
  f = g.fa
  K = get_order(f)
  np = length(x)
  ndofs = length(f)
  xi = testitem(x)
  Vr = _return_val_eltype(f,x)
  G = gradient_type(Vr,xi)
  a = zeros(G,(np,ndofs))
  mb = CartProdPolyBasis(PT,Val(D),V,K-1,_p_filter)
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
    # terms for ∇((𝓟ₖ)³)
    for j in 1:ndofsP
      a[i,j] = Px[i,j]
    end
    # terms for  ∇(x × (𝓟ₖ\𝓟ₖ₋₁)³)
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
    # terms for  ∇((𝓟ₖ)²)
    for j in 1:ndofsP
      a[i,j] = Px[i,j]
    end
    # terms for  ∇(x × (𝓟ₖ\𝓟ₖ₋₁)²)
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
