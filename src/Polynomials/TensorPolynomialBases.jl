#################################
# Tensorial nD polynomial bases #
#################################

"""
    struct TensorPolynomialBasis{D,V,K,PT} <: PolynomialBasis{D,V,K,PT}

Type representing a tensorial basis of (an)isotropic `D`-multivariate `V`-valued
(scalar, vector, or tensor) polynomial basis, that is:

The polynomial space is the tensor product of a scalar polynomial space (one for
each independant component of V).

The scalar polynomial basis is

  { x ⟶ b`ᴷ`\\_α(x) = b`ᴷ`\\_α₁(x₁) × b`ᴷ`\\_α₂(x₂) × ... × b`ᴷ`\\_α`D`(x`D`) |  α ∈ `terms` }

where b`ᴷ`\\_αᵢ(xᵢ) is the αᵢth 1D basis polynomial of the basis `PT` of order `K`
evaluated at xᵢ, and where α = (α₁, α₂, ..., α`D`) is a multi-index in `terms`,
a subset of ⟦0,`K`⟧`ᴰ`. `terms` is a field that can be passed in a constructor.

The fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to second
order derivatives.
"""
struct TensorPolynomialBasis{D,V,K,PT} <: PolynomialBasis{D,V,K,PT}
  orders::NTuple{D,Int}
  terms::Vector{CartesianIndex{D}}

  function TensorPolynomialBasis{D}(
    ::Type{PT},
    ::Type{V},
    orders::NTuple{D,Int},
    terms::Vector{CartesianIndex{D}}) where {D,V,PT<:Polynomial}

    K = maximum(orders; init=0)
    msg =  "Some term contain a higher index than the maximum degree + 1."
    @check all( term -> (maximum(Tuple(term), init=0) <= K+1), terms) msg
    new{D,V,K,PT}(orders,terms)
  end
end

@inline Base.size(a::TensorPolynomialBasis{D,V}) where {D,V} = (length(a.terms)*num_indep_components(V),)

function TensorPolynomialBasis(
   ::Type{PT},
   ::Val{D},
   ::Type{V},
   orders::NTuple{D,Int},
   terms::Vector{CartesianIndex{D}}) where {PT<:Polynomial,D,V}

  TensorPolynomialBasis{D}(PT,V,orders,terms)
end

"""
    TensorPolynomialBasis(::Type{PT}, ::Val{D}, ::Type{V}, orders::Tuple [, filter::Function])

This version of the constructor allows to pass a tuple `orders` containing the
polynomial order to be used in each of the `D` dimensions in order to  construct
an anisotropic tensor-product space.
"""
function TensorPolynomialBasis(
  ::Type{PT}, ::Val{D}, ::Type{V}, orders::NTuple{D,Int}, filter::Function=_q_filter
  ) where {PT,D,V}

  terms = _define_terms(filter, orders)
  TensorPolynomialBasis{D}(PT,V,orders,terms)
end

"""
    TensorPolynomialBasis(::Type{V}, ::Val{D}, order::Int [, filter::Function]) where {D,V}

Returns an instance of `TensorPolynomialBasis` representing a multivariate polynomial basis
in `D` dimensions, of polynomial degree `order`, whose value is represented by the type `V`.
The type `V` is typically `<:Number`, e.g., `Float64` for scalar-valued functions and `VectorValue{D,Float64}`
for vector-valued ones.

# Filter function

The `filter` function is used to select which terms of the tensor product space
of order `order` in `D` dimensions are to be used. If the filter is not provided, the full tensor-product
space is used by default leading to a multivariate polynomial space of type Q.
The signature of the filter function is

    (e,order) -> Bool

where `e` is a tuple of `D` integers containing the exponents of a multivariate monomial. The following filters
are used to select well known polynomial spaces

- Q space: `(e,order) -> true`
- P space: `(e,order) -> sum(e) <= order`
- "Serendipity" space: `(e,order) -> sum( [ i for i in e if i>1 ] ) <= order`

"""
function TensorPolynomialBasis(
  ::Type{PT}, VD::Val{D}, ::Type{V}, order::Int, filter::Function=_q_filter) where {PT,D,V}

  orders = tfill(order,VD)
  TensorPolynomialBasis(PT,Val(D),V,orders,filter)
end

# API

"""
    get_exponents(b::TensorPolynomialBasis)

Get a vector of tuples with the exponents of all the terms in the basis.

# Examples

```jldoctest
using Gridap.Polynomials

b = MonomialBasis(Val(2),Float64,2)

exponents = get_exponents(b)

println(exponents)

# output
Tuple{Int,Int}[(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2)]
```
"""
function get_exponents(b::TensorPolynomialBasis)
  indexbase = 1
  [Tuple(t) .- indexbase for t in b.terms]
end

"""
    get_orders(b::TensorPolynomialBasis)

Return the D-tuple of polynomial orders in each dimension
"""
function get_orders(b::TensorPolynomialBasis)
  b.orders
end


###########
# Helpers #
###########

_q_filter( e,order) = (maximum(e,init=0) <= order) # ℚₙ
_qs_filter(e,order) = (maximum(e,init=0) == order) # ℚₙ\ℚ₍ₙ₋₁₎
_p_filter( e,order) = (sum(e) <= order)     # ℙₙ
_ps_filter(e,order) = (sum(e) == order)     # ℙₙ\ℙ₍ₙ₋₁₎

function _p_dim(order,D)
  dim = 1
  for d in 1:D
    dim *= order+d
  end
  dim/factorial(D)
end

function _define_terms(filter,orders)
  t = orders .+ 1
  g = (0 .* orders) .+ 1
  cis = CartesianIndices(t)
  co = CartesianIndex(g)
  maxorder = maximum(orders, init=0)
  [ ci for ci in cis if filter(Int[Tuple(ci-co)...],maxorder) ]
end


#################################
# nD evaluations implementation #
#################################

function _evaluate_nd!(
  b::TensorPolynomialBasis{D,V,K,PT}, x,
  r::AbstractMatrix{V}, i,
  v::AbstractVector{V},
  c::AbstractMatrix{T}) where {D,V,K,PT,T}

  terms  = b.terms
  orders = b.orders

  for d in 1:D
    Kd = Val(orders[d])
    _evaluate_1d!(PT,Kd,c,x,d)
  end

  k = 1
  for ci in terms

    s = one(T)
    for d in 1:D
      @inbounds s *= c[d,ci[d]]
    end

    k = _tensorial_set_value!(v,s,k)
  end

  #r[i] = v
  @inbounds for j in 1:length(b)
      r[i,j] = v[j]
  end
end

"""
    _tensorial_set_value!(v::AbstractVector{<:Real},s,k)

v[k] = s; return k+1
"""
function _tensorial_set_value!(v::AbstractVector{<:Real},s,k)
    @inbounds v[k] = s
    k+1
end

"""
    _tensorial_set_value!(v::AbstractVector{V},s::T,k)

```
v[k]   = V(s, 0, ..., 0)
v[k+1] = V(0, s, ..., 0)
        ⋮
v[k+N] = V(0, ..., 0, s)
return k+N
```

where `N` is the number of independent components of `V<:MultiValue`.

This means that the basis has the same polynomial space in each component, so it
is tensorial relative to the `V` components.
"""
function _tensorial_set_value!(v::AbstractVector{V},s::T,k) where {V,T}
  ncomp = num_indep_components(V)
  z = zero(T)
  @inbounds for j in 1:ncomp
    v[k] = ntuple(i -> ifelse(i == j, s, z),Val(ncomp))
    k += 1
  end
  k
end

function _gradient_nd!(
  b::TensorPolynomialBasis{D,V,K,PT}, x,
  r::AbstractMatrix{G}, i,
  v::AbstractVector{G},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  s::MVector{D,T}) where {D,V,K,PT,G,T}

  terms  = b.terms
  orders = b.orders

  for d in 1:D
    Kd = Val(orders[d])
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

    k = _tensorial_set_gradient!(v,s,k,V)
  end

  #r[i] = v
  @inbounds for j in 1:length(b)
      r[i,j] = v[j]
  end
end

function _tensorial_set_gradient!(
  v::AbstractVector{G},s,k,::Type{<:Real}) where G

  @inbounds v[k] = s
  k+1
end

# TODO comment
@generated function _tensorial_set_gradient!(
  v::AbstractVector{G},s,k,::Type{V}) where {G,V}
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
    body *= "@inbounds v[k] = ($(join(tuple(m...), ", ")));"
    body *= "k = k + 1;"
  end

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k))
end

# Specialization for SymTensorValue and SymTracelessTensorValue,
# necessary as long as outer(Point, V<:AbstractSymTensorValue)::G does not
# return a tensor type G that implements the appropriate symmetries of the
# gradient (and hessian)
#
# This is still (independant-)component tensorial as each independent SymTensor
# component holds the same (scalar multivariate) polynomial space.
@generated function _tensorial_set_gradient!(
  v::AbstractVector{G},s,k,::Type{V}) where {G,V<:AbstractSymTensorValue{D}} where D
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
      body *= "@inbounds v[k] = ($(join(tuple(m...), ", ")));"
      body *= "k = k + 1;"
    end
  end

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k))
end

function _hessian_nd!(
  b::TensorPolynomialBasis{D,V,K,PT}, x,
  r::AbstractMatrix{G}, i,
  v::AbstractVector{G},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  s::MMatrix{D,D,T}) where {D,V,K,PT,G,T}

  terms  = b.terms
  orders = b.orders

  for d in 1:D
    Kd = Val(orders[d])
    _derivatives_1d!(PT,Kd,(c,g,h),x,d)
  end

  k = 1

  for ci in terms

    for i in eachindex(s)
      @inbounds s[i] = one(T)
    end

    for r in 1:D
      for q in 1:D
        for d in 1:D
          if d != q && d != r
            @inbounds s[r,q] *= c[d,ci[d]]
          elseif d == q && d ==r
            @inbounds s[r,q] *= h[d,ci[d]]
          else
            @inbounds s[r,q] *= g[d,ci[d]]
          end
        end
      end
    end

    k = _tensorial_set_gradient!(v,s,k,V)
  end

  #r[i] = v
  @inbounds for j in 1:length(b)
      r[i,j] = v[j]
  end
end

