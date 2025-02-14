#################################
# Tensorial nD polynomial bases #
#################################

"""
    struct CartProdPolyBasis{D,V,K,PT} <: PolynomialBasis{D,V,K,PT}

"Cartesian product polynomial basis"

Type representing a basis of a (an)isotropic `D`-multivariate `V`-valued
cartesian product polynomial space

`V`(𝕊, 𝕊, ..., 𝕊)

where 𝕊 is a scalar multivariate polynomial space.

The scalar polynomial basis spanning 𝕊 is defined as

  { x ⟶ bα`ᴷ`(x) = bα₁`ᴷ`(x₁) × bα₂`ᴷ`(x₂) × ... × bα`Dᴷ`(x`D`) |  α ∈ `terms` }

where bαᵢ`ᴷ`(xᵢ) is the αᵢth 1D basis polynomial of the basis `PT` of order `K`
evaluated at xᵢ (iᵗʰ comp. of x), and where α = (α₁, α₂, ..., α`D`) is a
multi-index in `terms`, a subset of ⟦0,`K`⟧`ᴰ`. `terms` is a field that can be
passed in a constructor.

This type fully implements the [`Field`](@ref) interface, with up to second
order derivatives.
"""
struct CartProdPolyBasis{D,V,K,PT} <: PolynomialBasis{D,V,K,PT}
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
    new{D,V,K,PT}(orders,terms)
  end
end

@inline Base.size(a::CartProdPolyBasis{D,V}) where {D,V} = (length(a.terms)*num_indep_components(V),)

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
  Tuple(Tuple(t) .- indexbase for t in b.terms)
end

"""
    get_orders(b::CartProdPolyBasis)

Return the D-tuple of polynomial orders in each spatial dimension
"""
function get_orders(b::CartProdPolyBasis)
  b.orders
end

#################################
# nD evaluations implementation #
#################################

function _evaluate_nd!(
  b::CartProdPolyBasis{D,V,K,PT}, x,
  r::AbstractMatrix{V}, i,
  c::AbstractMatrix{T}) where {D,V,K,PT,T}

  terms  = b.terms
  orders = b.orders

  Kd = Val(K)
  for d in 1:D
    # The optimization below of fine tuning Kd is a bottlneck if not put in a
    # function due to runtime dispatch and creation of Val(Kd)
    # Kd = Val(orders[d])
    _evaluate_1d!(PT,Kd,c,x,d)
  end

  k = 1
  for ci in terms

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

where `N = num_indep_components(V)`.
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
  b::CartProdPolyBasis{D,V,K,PT}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  s::MVector{D,T}) where {D,V,K,PT,G,T}

  terms  = b.terms
  orders = b.orders

  Kd = Val(K)
  for d in 1:D
    _derivatives_1d!(PT,Kd,(c,g),x,d)
  end

  k = 1
  for ci in terms

    for i in eachindex(s)
      @inbounds s[i] = one(T)
    end

    for q in 1:D
      for d in 1:D
        if d != q
          @inbounds s[q] *= c[d,ci[d]]
        else
          @inbounds s[q] *= g[d,ci[d]]
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
"""
@generated function _cartprod_set_derivative!(
  r::AbstractMatrix{G},i,s,k,::Type{V}) where {G,V}
  # Git blame me for readable non-generated version

  w = zero(V)
  m = Array{String}(undef, size(G))
  N_val_dims = length(size(V))
  s_size = size(G)[1:end-N_val_dims]

  body = "T = eltype(s); z = zero(T);"
  for ci in CartesianIndices(s_size)
    id = join(Tuple(ci))
    body *= "@inbounds s$id = s[$ci];"
  end

  for j in CartesianIndices(w)
    for i in CartesianIndices(m)
      m[i] = "z"
    end
    for ci in CartesianIndices(s_size)
      id = join(Tuple(ci))
      m[ci,j] = "s$id"
    end
    body *= "@inbounds r[i,k] = ($(join(tuple(m...), ", ")));"
    body *= "k = k + 1;"
  end

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k))
end

# Specialization for SymTensorValue and SymTracelessTensorValue,
# necessary as long as outer(Point, V<:AbstractSymTensorValue)::G does not
# return a tensor type G that implements the appropriate symmetries of the
# gradient (and hessian)
@generated function _cartprod_set_derivative!(
  r::AbstractMatrix{G},i,s,k,::Type{V}) where {G,V<:AbstractSymTensorValue{D}} where D
  # Git blame me for readable non-generated version

  T = eltype(s)
  m = Array{String}(undef, size(G))
  s_length = size(G)[1]

  is_traceless = V <: SymTracelessTensorValue
  skip_last_diagval = is_traceless ? 1 : 0    # Skid V_DD if traceless

  body = "z = $(zero(T));"
  for i in 1:s_length
    body *= "@inbounds s$i = s[$i];"
  end

  for c in 1:(D-skip_last_diagval) # Go over cols
    for r in c:D                   # Go over lower triangle, current col
      for i in eachindex(m)
        m[i] = "z"
      end
      for i in 1:s_length # indices of the Vector s
        m[i,r,c] = "s$i"
        if (r!=c)
          m[i,c,r] = "s$i"
        elseif is_traceless # V_rr contributes negatively to V_DD (tracelessness)
          m[i,D,D] = "-s$i"
        end
      end
      body *= "@inbounds r[i,k] = ($(join(tuple(m...), ", ")));"
      body *= "k = k + 1;"
    end
  end

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k))
end

function _hessian_nd!(
  b::CartProdPolyBasis{D,V,K,PT}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  s::MMatrix{D,D,T}) where {D,V,K,PT,G,T}

  terms  = b.terms
  orders = b.orders

  Kd = Val(K)
  for d in 1:D
    _derivatives_1d!(PT,Kd,(c,g,h),x,d)
  end

  k = 1

  for ci in terms

    for i in eachindex(s)
      @inbounds s[i] = one(T)
    end

    for t in 1:D
      for q in 1:D
        for d in 1:D
          if d != q && d != t
            @inbounds s[t,q] *= c[d,ci[d]]
          elseif d == q && d ==t
            @inbounds s[t,q] *= h[d,ci[d]]
          else
            @inbounds s[t,q] *= g[d,ci[d]]
          end
        end
      end
    end

    k = _cartprod_set_derivative!(r,i,s,k,V)
  end
end

