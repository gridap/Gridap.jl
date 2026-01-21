#################################
# Tensorial nD polynomial bases #
#################################

"""
    struct CartProdPolyBasis{D,V,PT} <: PolynomialBasis{D,V,PT}

"Cartesian product of a scalar tensor polynomial basis"

Type representing a basis of a (an)isotropic `D`-multivariate `V`-valued
cartesian product polynomial space

`V`(𝕊, ∅, ..., ∅) ⊕ `V`(∅, 𝕊, ∅, ..., ∅) ⊕ ... ⊕ `V`(∅, ..., ∅, 𝕊)

where the scalar space 𝕊 is a (subspace of a) tensor product space of an
univariate polynomial basis.

The scalar polynomial basis spanning 𝕊 is defined as

  { x ⟶ bα`ᴷ`(x) = bα₁`ᴷ`(x₁) × bα₂`ᴷ`(x₂) × ... × bα`Dᴷ`(x`D`) |  α ∈ `terms` }

where bαᵢ`ᴷ`(xᵢ) is the αᵢth 1D basis polynomial of the basis `PT` of order `K`
evaluated at xᵢ (iᵗʰ comp. of x), and where α = (α₁, α₂, ..., α`D`) is a
multi-index in `terms`, a subset of {0:`K`}`ᴰ`. `terms` is a field that can be
passed in a constructor.

This type fully implements the [`Field`](@ref) interface, with up to second
order derivatives.
"""
struct CartProdPolyBasis{D,V,PT} <: PolynomialBasis{D,V,PT}
  max_order::Int
  orders::NTuple{D,Int}
  terms::Vector{CartesianIndex{D}}

  function CartProdPolyBasis{D}(
    ::Type{PT},
    ::Type{V},
    orders::NTuple{D,Int},
    terms::Vector{CartesianIndex{D}}) where {D,V,PT<:Polynomial}

    @check isconcretetype(PT) "PT needs to be a concrete <:Polynomial type"

    K = maximum(orders; init=0)
    msg =  "Some term contain a higher index than the maximum degree + 1."
    @check all( term -> (maximum(Tuple(term), init=0) <= K+1), terms) msg
    new{D,V,PT}(K,orders,terms)
  end
end

@inline Base.size(a::CartProdPolyBasis{D,V}) where {D,V} = (length(a.terms)*num_indep_components(V),)
@inline get_order(b::CartProdPolyBasis) = b.max_order

function testvalue(::Type{CartProdPolyBasis{D,V,PT}}) where {D,V,PT}
  CartProdPolyBasis{D}(PT,V,tfill(0,Val(D)),CartesianIndex{D}[])
end

function CartProdPolyBasis(
   ::Type{PT},
   ::Val{D},
   ::Type{V},
   orders::NTuple{D,Int},
   terms::Vector{CartesianIndex{D}}) where {PT<:Polynomial,D,V}

  CartProdPolyBasis{D}(PT,V,orders,terms)
end

"""
    CartProdPolyBasis(::Type{PT}, ::Val{D}, ::Type{V}, orders::Tuple [, filter=_q_filter])

This constructor allows to pass a tuple `orders` containing the maximal
polynomial order to be used in each of the `D` spatial dimensions in order to
construct a tensorial anisotropic `D`-multivariate space 𝕊.

If a filter is provided, it is applied on the cartesian product terms
CartesianIndices(`orders`), with maximum(`orders`) as order argument.
"""
function CartProdPolyBasis(
  ::Type{PT}, ::Val{D}, ::Type{V}, orders::NTuple{D,Int}, filter::Function=_q_filter
  ) where {PT,D,V}

  terms = _define_terms(filter, orders)
  CartProdPolyBasis{D}(PT,V,orders,terms)
end

"""
    CartProdPolyBasis(::Type{PT}, ::Type{V}, ::Val{D}, order::Int [, filter=_q_filter])

Return a `CartProdPolyBasis{D,V,order,PT}` where 𝕊 is defined by the terms
filtered by

    term -> filter(term, order).

See the [Filter functions](@ref) section of the documentation for more details.
"""
function CartProdPolyBasis(
  ::Type{PT}, VD::Val{D}, ::Type{V}, order::Int, filter::Function=_q_filter) where {PT,D,V}

  orders = tfill(order,VD)
  CartProdPolyBasis(PT,Val(D),V,orders,filter)
end

# API

"""
    get_exponents(b::CartProdPolyBasis)

Get a vector of tuples with the exponents of all the terms in the basis of 𝕊,
the components scalar space of `b`.

# Example

```jldoctest
using Gridap.Polynomials

b = MonomialBasis(Val(2),Float64,2)

exponents = get_exponents(b)

println(exponents)

# output
Tuple{Int,Int}[(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2)]
```
"""
function get_exponents(b::CartProdPolyBasis)
  indexbase = 1
  [ Tuple(t) .- indexbase for t in b.terms ]
end

function get_orders(b::CartProdPolyBasis)
  b.orders
end

#################################
# nD evaluations implementation #
#################################

function _evaluate_nd!(
  b::CartProdPolyBasis{D,V,PT}, x,
  r::AbstractMatrix, i,
  c::AbstractMatrix{T}, K) where {D,V,PT,T}

  for d in 1:D
    Kd = b.orders[d]
    _evaluate_1d!(PT,Kd,c,x,d)
  end

  k = 1
  for ci in b.terms

    s = one(T)
    for d in 1:D
      @inbounds s *= c[d,ci[d]]
    end

    k = _cartprod_set_value!(r,i,s,k)
  end
end

"""
    _cartprod_set_value!(r::AbstractMatrix{<:Real},i,s,k)

    r[i,k] = s; return k+1

`s` is scalar
"""
function _cartprod_set_value!(r::AbstractMatrix{<:Real},i,s,k)
  @inbounds r[i,k] = s
  k+1
end

"""
    _cartprod_set_value!(r::AbstractMatrix{V},i,s::T,k,l)

```
r[i,k]     = V(s, 0,    ..., 0)
r[i,k+1]   = V(0, s, 0, ..., 0)
⋮
r[i,k+N-1] = V(0,    ..., 0, s)
return k+N
```

where `N = num_indep_components(V)`, and `s` is scalar.
"""
function _cartprod_set_value!(r::AbstractMatrix{V},i,s::T,k) where {V,T}
  ncomp = num_indep_components(V)
  z = zero(T)
  @inbounds for j in 1:ncomp
    r[i,k] = ntuple(i -> ifelse(i == j, s, z),Val(ncomp))
    k += 1
  end
  k
end

function _gradient_nd!(
  b::CartProdPolyBasis{D,V,PT}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  s::MVector{D,T}, K) where {D,V,PT,G,T}

  for d in 1:D
    _derivatives_1d!(PT,K,(c,g),x,d)
  end

  k = 1
  @inbounds for ci in b.terms

    s[:] .= one(T)

    for q in 1:D
      for d in 1:D
        if d != q
          s[q] *= c[d,ci[d]]
        else
          s[q] *= g[d,ci[d]]
        end
      end
    end

    k = _cartprod_set_derivative!(r,i,s,k,V)
  end
end

"""
    _cartprod_set_derivative!(r::AbstractMatrix{G},i,s,k,::Type{<:Real})

```
r[i,k] = s = (∇bᵏ)(xi); return k+1
```

where bᵏ is the kᵗʰ basis polynomial. Note that `r[i,k]` is a `VectorValue` or
`TensorValue` and `s` a `MVector` or `MMatrix` respectively, of same size.
"""
function _cartprod_set_derivative!(
  r::AbstractMatrix{G},i,s,k,::Type{<:Real}) where G

  @inbounds r[i,k] = s
  k+1
end

"""
    _cartprod_set_derivative!(r::AbstractMatrix{G},i,s,k,::Type{V})

```
z = zero(s)
r[i,k]     = G(s…, z…,     ..., z…) = (Dbᵏ    )(xi)
r[i,k+1]   = G(z…, s…, z…, ..., z…) = (Dbᵏ⁺¹  )(xi)
⋮
r[i,k+n-1] = G(z…,     ..., z…, s…) = (Dbᵏ⁺ⁿ⁻¹)(xi)
return k+n
```

Note that `r[i,k]` is a `TensorValue` or `ThirdOrderTensorValue` and `s` a
`MVector` or `MMatrix`.

# Extended help

When `V` has dependent components, the derivative values are not defined as
simply as above, but this function implements what it should: seting derivatives
with respect to each independent components of `V`.
"""
@generated function _cartprod_set_derivative!(
  r::AbstractMatrix{G},i,s,k,::Type{V}) where {G,V<:MultiValue}

  T = eltype(s)
  zs = zero(s)
  m = Array{String}(undef, size(G))

  V_deri = MultiValue(zs)
  grad_V_basis = [ ∇i⊗vj for ∇i in component_basis(V_deri), vj in component_basis(V)]

  body = "z = $(zero(T));"
  for i in 1:length(zs)
    body *= "@inbounds s$i = s[$i];"
  end

  for icomp in 1:num_indep_components(V)
    m .= "z"
    for i in 1:length(zs)
      gVi_comps = grad_V_basis[i,icomp]
      for j in 1:length(G)
        gVij = gVi_comps[j]
        if !iszero(gVij)
          #m[j] += grad_V_basis[i,icomp][j] * s[i]
          m[j] *= " + $gVij*s$i"
        end
      end
    end
    # r[i,k] = sum( grad_V_basis[i,icomp] .* s[i] for i in 1:length(s) )
    body *= "@inbounds r[i,k] = ($(join(tuple(m...), ", ")));"
    body *= "k = k + 1;"
  end

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k))
end

function _hessian_nd!(
  b::CartProdPolyBasis{D,V,PT}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  s::MMatrix{D,D,T}, K) where {D,V,PT,G,T}

  for d in 1:D
    _derivatives_1d!(PT,K,(c,g,h),x,d)
  end

  k = 1

  @inbounds for ci in b.terms

    s[:] = one(T)

    for t in 1:D
      for q in 1:D
        for d in 1:D
          if d != q && d != t
            s[t,q] *= c[d,ci[d]]
          elseif d == q && d ==t
            s[t,q] *= h[d,ci[d]]
          else
            s[t,q] *= g[d,ci[d]]
          end
        end
      end
    end

    k = _cartprod_set_derivative!(r,i,s,k,V)
  end
end

