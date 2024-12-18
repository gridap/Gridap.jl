"""
    CompWiseTensorPolyBasis{D,V,K,PT,L} <: PolynomialBasis{D,V,K,PT}

"Polynomial basis of component wise tensor product spaces"

Polynomial basis for a multivariate `MultiValue`'d polynomial space
`V`(𝕊¹, 𝕊², ..., 𝕊ᴸ) with `L`>1, where the scalar multivariate spaces 𝕊ˡ
(for 1 ≤ l ≤ `L`) of each (independent) component of `V` is the tensor product
of 1D ℙ¹ spaces of order oₗₙ for 1 ≤ n ≤ `D`, that is:

𝕊¹ = ℙ¹ₒ(1,1) ⊗ … ⊗ ℙ¹ₒ(1,`D`)\\
⋮\\
𝕊ˡ =     ⊗ₙ  ℙ¹ₒ₍ₗ,ₙ₎\\
⋮\\
𝕊ᴸ = ℙ¹ₒ(`L`,1) ⊗ … ⊗ ℙ¹ₒ(`L`,`D`)

The `L`×`D` matrix of orders o is given in the constructor, and `K` is the
maximum of o.
"""
struct CompWiseTensorPolyBasis{D,V,K,PT,L} <: PolynomialBasis{D,V,K,PT}
  orders::SMatrix{L,D,Int}

  function CompWiseTensorPolyBasis{D}(
    ::Type{PT}, ::Type{V}, orders::SMatrix{L,D,Int}) where {D,PT<:Polynomial,V,L}

    msg1 = "The orders matrix rows number must match the number of independent components of V"
    @check L == num_indep_components(V) msg1
    msg2 = "The Component Wise construction is useless for one component, use TensorPolynomialBasis instead"
    @check L > 1 msg2
    @check D > 0
    K = maximum(orders)

    new{D,V,K,PT,L}(orders)
  end
end

Base.size(a::CompWiseTensorPolyBasis) = ( sum(prod.(eachrow(a.orders .+ 1))), )

"""
    get_comp_terms(f::CompWiseTensorPolyBasis)

Return a `NTuple{L,CartesianIndices{D}}` containing, for each component
1 ≤ l ≤ `L`, the Cartesian indices iterator over the terms
in ⟦1,`o`(l,1)+1⟧ × ⟦1,`o`(l,2)+1⟧ × … × ⟦1,`o`(l,D)+1⟧ that define 𝕊ˡ.

E.g., if `orders=[ 0 1; 1 0]`, then the `comp_terms` are
`( CartesianIndices{2}((1,2)), CartesianIndices{2}((2,1)) )`.
"""
function get_comp_terms(f::CompWiseTensorPolyBasis{D,V,K,PT,L}) where {D,V,K,PT,L}
  _terms(l) = CartesianIndices( Tuple(f.orders[l,:] .+ 1) )
  comp_terms = ntuple(l -> _terms(l), Val(L))
  comp_terms::NTuple{L,CartesianIndices{D}}
end


########################
# Field implementation #
########################

function evaluate!(
  cache,
  f::CompWiseTensorPolyBasis{D,V,K,PT,L},
  x::AbstractVector{<:Point}) where {D,V,K,PT,L}

  r, v, c = cache
  np = length(x)
  ndof = length(f)
  ndof_1d = get_order(f) + 1 # K + 1
  comp_terms = get_comp_terms(f)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,ndof_1d))
  for i in 1:np
    @inbounds xi = x[i]
    _evaluate_nd_cwtpb!(PT,v,xi,f.orders,comp_terms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function evaluate!(
  cache,
  fg::FieldGradientArray{1,<:CompWiseTensorPolyBasis{D,V,K,PT}},
  x::AbstractVector{<:Point}) where {D,V,K,PT}

  f = fg.fa
  r, v, c, g = cache
  np = length(x)
  ndof = length(f)
  ndof_1d = get_order(f) + 1
  comp_terms = get_comp_terms(f)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,ndof_1d))
  setsize!(g,(D,ndof_1d))
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_cwtpb!(PT,v,xi,f.orders,comp_terms,c,g,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function evaluate!(
  cache,
  fg::FieldGradientArray{2,<:CompWiseTensorPolyBasis{D,V,K,PT}},
  x::AbstractVector{<:Point}) where {D,V,K,PT}

  f = fg.fa
  r, v, c, g, h = cache
  np = length(x)
  ndof = length(f)
  ndof_1d = get_order(f) + 1
  comp_terms = get_comp_terms(f)
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,ndof_1d))
  setsize!(g,(D,ndof_1d))
  setsize!(h,(D,ndof_1d))
  for i in 1:np
    @inbounds xi = x[i]
    _hessian_nd_cwtpb!(PT,v,xi,f.orders,comp_terms,c,g,h,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end


###########
# Helpers #
###########

function _evaluate_nd_cwtpb!(
  PT::Type{<:Polynomial},
  v::AbstractVector{V},
  x,
  orders::SMatrix{L,D,Int},
  comp_terms::NTuple{L,CartesianIndices{D}},
  c::AbstractMatrix{T}) where {V,L,D,T}

  for d in 1:D
    # for each coordinate d, the order at which the basis should be evaluated is
    # the maximum d-order for any component l
    Kd = Val(maximum(orders[:,d]))
    _evaluate_1d!(PT,Kd,c,x,d)
  end

  m = zero(Mutable(V))
  k = 1

  for (l,terms) in enumerate(comp_terms)
    for ci in terms

      s = one(T)
      @inbounds for d in 1:D
        s *= c[d,ci[d]]
      end

      k = _comp_wize_set_value!(v,s,k,l)
    end
  end
end

"""
    _comp_wize_set_value!(v::AbstractVector{V},s::T,k,l)

```
v[k]   = V(0, ..., 0, s, 0, ..., 0); return k+1
```

where `s` is at position `l` in `V<:MultiValue`.
"""
function _comp_wize_set_value!(v::AbstractVector{V},s::T,k,l) where {V,T}
  z = zero(T)
  ncomp = num_indep_components(V)
  v[k] = ntuple(i -> ifelse(i == l, s, z),Val(ncomp))
  return k + 1
end

function _gradient_nd_cwtpb!(
  PT::Type{<:Polynomial},
  v::AbstractVector{G},
  x,
  orders::SMatrix{L,D,Int},
  comp_terms::NTuple{L,CartesianIndices{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {G,L,D,T,V}

  for d in 1:D
    # for each spatial coordinate d, the order at which the basis should be
    # evaluated is the maximum d-order for any component l
    Kd = Val(maximum(orders[:,d]))
    _derivatives_1d!(PT,Kd,(c,g),x,d)
  end

  s = zero(Mutable(V))
  k = 1

  for (l,terms) in enumerate(comp_terms)
    for ci in terms

      for i in eachindex(s)
        s[i] = one(T)
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

      k = _comp_wize_set_gradient!(v,s,k,Val(l),V)
    end
  end
end

function _comp_wize_set_gradient!(
  v::AbstractVector{G},s,k,l,::Type{<:Real}) where G

  @inbounds v[k] = s
  k+1
end

@generated function _comp_wize_set_gradient!(
  v::AbstractVector{G},s,k,::Val{l},::Type{V}) where {G,l,V}

  N_val_dims = length(size(V))
  s_size = size(G)[1:end-N_val_dims]

  body = "T = eltype(s); z = zero(T);"
  m = Array{String}(undef, size(G))
  m .= "z"

  for ci in CartesianIndices(s_size)
    m[ci,l] = "(@inbounds s[$ci])"
  end
  body *= "@inbounds v[k] = ($(join(tuple(m...), ", ")));"

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k+1))
end

@generated function _comp_wize_set_gradient!(
  v::AbstractVector{G},s,k,::Type{V}) where {G,V<:AbstractSymTensorValue{D}} where D

  @notimplemented
end

function _hessian_nd_cwtpb!(
  PT::Type{<:Polynomial},
  v::AbstractVector{G},
  x,
  orders::SMatrix{L,D,Int},
  comp_terms::NTuple{L,CartesianIndices{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  ::Type{V}) where {G,L,D,T,V}

  for d in 1:D
    # for each spatial coordinate d, the order at which the basis should be
    # evaluated is the maximum d-order for any component l
    Kd = Val(maximum(orders[:,d]))
    _derivatives_1d!(PT,Kd,(c,g,h),x,d)
  end

  s = zero(Mutable(V))
  k = 1

  for (l,terms) in enumerate(comp_terms)
    for ci in terms

      for i in eachindex(s)
        s[i] = one(T)
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

      k = _comp_wize_set_gradient!(v,s,k,l,V)
    end
  end
end


################################
# Basis for Nedelec on D-cubes #
################################

"""
    QGradBasis(::Type{PT}, ::Val{D}, ::Type{T}, order::Int) :: PolynomialBasis

Return a basis of (ℚ\\_order)ᴰ ⊕ x × ( ℚ\\_order \\ ℚ\\_{order-1}) )ᴰ, the
polynomial space for Nedelec elements on `D`-dimensonal cubes with scalar type `T`.

The `order` argument has the following meaning: the curl of the  functions in
this basis is in the ℚ space of degree `order`.

`PT<:Polynomial` is the choice of scalar 1D polynomial basis.

For more details, see [`CompWiseTensorPolyBasis`](@ref), as `QGradBasis` returns
an instance of `CompWiseTensorPolyBasis{D,VectorValue{D,T},order+1,PT,D}`.
"""
function QGradBasis(::Type{PT},::Val{D},::Type{T},order::Int) where {PT,D,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"

  V = VectorValue{D,T}
  m = [ order + (i==j ? 0 : 1) for i in 1:D, j in 1:D ]
  orders = SMatrix{D,D,Int}(m)
  CompWiseTensorPolyBasis{D}(PT, V, orders)
end

function QGradBasis(::Type{PT},::Val{1},::Type{T},order::Int) where {PT,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"

  V = VectorValue{1,T}
  TensorPolynomialBasis(PT, Val(1), V, order+1)
end


#######################################
# Basis for Raviart-Thomas on D-cubes #
#######################################

"""
    QCurlGradBasis(::Type{PT}, ::Val{D}, ::Type{T}, order::Int) :: PolynomialBasis

Return a basis of (ℚ\\_order)ᴰ ⊕ x (ℚ\\_order \\ ℚ\\_{order-1}), the polynomial space
for Raviart-Thomas elements on `D`-dimensonal cubes with scalar type `T`.

The `order` argument has the following meaning: the divergence of the functions
in this basis is in the ℚ space of degree `order`.

`PT<:Polynomial` is the choice of scalar 1D polynomial basis.

For more details, see [`CompWiseTensorPolyBasis`](@ref), as `QCurlGradBasis` returns
an instance of `CompWiseTensorPolyBasis{D,VectorValue{D,T},order+1,PT,D}`.
"""
function QCurlGradBasis(::Type{PT},::Val{D},::Type{T},order::Int) where {PT,D,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"

  V = VectorValue{D,T}
  m = [ order + (i==j ? 1 : 0) for i in 1:D, j in 1:D ]
  orders = SMatrix{D,D,Int}(m)
  CompWiseTensorPolyBasis{D}(PT, V, orders)
end

function QCurlGradBasis(::Type{PT},::Val{1},::Type{T},order::Int) where {PT,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"

  V = VectorValue{1,T}
  TensorPolynomialBasis(PT, Val(1), V, order+1)
end
