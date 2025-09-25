"""
    CompWiseTensorPolyBasis{D,V,PT} <: PolynomialBasis{D,V,PT}

"Polynomial basis of component wise tensor product polynomial spaces"

Polynomial basis for a `D`-multivariate `V`-valued polynomial space:

`V`(𝕊¹, 𝕊², ..., 𝕊ᴸ)

with `L`>1, where the scalar `D`-multivariate spaces 𝕊ˡ (for 1 ≤ l ≤ `L`) of each
(independent) component of `V` is the tensor product of 1D ℙ spaces of order
α(l,n) for 1 ≤ n ≤ `D`, that is:

𝕊¹ = ℙα(1,1) ⊗ … ⊗ ℙα(1,`D`)\\
⋮\\
𝕊ˡ =     ⊗ₙ  ℙα(l,n)\\
⋮\\
𝕊ᴸ = ℙα(`L`,1) ⊗ … ⊗ ℙα(`L`,`D`)

The `L`×`D` matrix of orders α is given in the constructor, and `K` is the
maximum of α. Any 1D polynomial family `PT<:Polynomial` is usable.

# Examples
These return instances of `CompWiseTensorPolyBasis`
```jldoctest
# a basis for Raviart-Thomas on quadrilateral with divergence in ℚ₃
b = FEEC_poly_basis(Val(2),Float64,4,1,:Q⁻; rotate_90)

# a basis for Raviart-Thomas on hexahedra with divergence in ℚ₃
b = FEEC_poly_basis(Val(3),Float64,4,2,:Q⁻)

# a basis for Nedelec on triangle with curl in ℚ₃
b = FEEC_poly_basis(Val(2),Float64,4,1,:Q⁻)

# a basis for Nedelec on hexahedra with curl in ℚ₃
b = FEEC_poly_basis(Val(3),Float64,4,1,:Q⁻)
```
"""
struct CompWiseTensorPolyBasis{D,V,PT} <: PolynomialBasis{D,V,PT}
  max_order::Int
  orders::Matrix{Int}
  comp_terms::Vector{CartesianIndices{D,NTuple{D,Base.OneTo{Int}}}} # of length L

  @doc"""
    CompWiseTensorPolyBasis{D}(PT,V,orders::Matrix{Int})

  `PT` is a [`Polynomial `](@ref) type, `V` a number type, and `orders` a L×D
  matrix where L is the number of independent components of `V`.
  """
  function CompWiseTensorPolyBasis{D}(
    ::Type{PT}, ::Type{V}, orders::Matrix{Int}) where {D,PT<:Polynomial,V}

    @check D > 0
    L = size(orders,1)
    msg1 = "The orders matrix rows number must match the number of independent components of V"
    @check L == num_indep_components(V) msg1
    msg2 = "The Component Wise construction is useless for one component, use CartProdPolyBasis instead"
    @check L > 1 msg2
    msg3 = "The orders matrix column number must match the number of spatial dimensions"
    @check size(orders,2) == D msg3
    @check isconcretetype(PT) "PT needs to be a concrete <:Polynomial type"

    K = maximum(orders)
    comp_terms = _compute_comp_terms(Val(D),V,orders)

    new{D,V,PT}(K,orders,comp_terms)
  end
end

Base.size(a::CompWiseTensorPolyBasis) = ( sum(prod.(eachrow(a.orders .+ 1))), )
get_order(b::CompWiseTensorPolyBasis) = b.max_order

function testvalue(::Type{<:CompWiseTensorPolyBasis{D,V,PT}}) where {D,V,PT}
  L = num_indep_components(V)
  orders = zero(Matrix{Int}(undef, (L,D)))
  CompWiseTensorPolyBasis{D}(PT,V,orders)
end

"""
    get_comp_terms(f::CompWiseTensorPolyBasis{D,V})

Return a tuple (terms\\_1, ..., terms\\_l, ..., terms\\_L) containing, for each
component of V, the Cartesian indices iterator over the terms that define 𝕊ˡ,
that is all elements of {1 : `o`(l,1)+1} × {1 : `o`(l,2)+1} × … × {1 : `o`(l,D)+1}.

E.g., if `orders=[ 0 1; 1 0]`, then the `comp_terms` are
`( CartesianIndices{2}((1,2)), CartesianIndices{2}((2,1)) )`.
"""
get_comp_terms(f::CompWiseTensorPolyBasis) = f.comp_terms

function _compute_comp_terms(::Val{D},::Type{V},orders) where {D,V}
  L = num_indep_components(V)
  _terms(l) = CartesianIndices( Tuple(orders[l,:] .+ 1)::NTuple{D,Int} )
  [ _terms(l) for l in 1:L ]
end


#################################
# nD evaluations implementation #
#################################

function _evaluate_nd!(
  b::CompWiseTensorPolyBasis{D,V,PT}, x,
  r::AbstractMatrix, i,
  c::AbstractMatrix{T}, ::Val{K}) where {D,V,PT,T,K}

  # optimization if PT is hierarchical: lower order polynomials do not depend on the maximum order
  if isHierarchical(PT)
    for d in 1:D
      _evaluate_1d!(PT,K,c,x,d)
    end
  end

  k = 1
  @inbounds for (l,terms) in enumerate(get_comp_terms(b))

    if !isHierarchical(PT)
      for d in 1:D
        # compute 1D polynomials for first component or recompute them if order changed
        kld = b.orders[l,d]
        if isone(l) || kld ≠ b.orders[l-1,d]
          _evaluate_1d!(PT,kld,c,x,d)
        end
      end
    end

    for ci in terms
      s = one(T)

      for d in 1:D
        s *= c[d,ci[d]]
      end

      k = _comp_wize_set_value!(r,i,s,k,l)
    end
  end
end

"""
    _comp_wize_set_value!(r::AbstractMatrix{V},i,s::T,k,l)

```
r[i,k]   = V(0, ..., 0, s, 0, ..., 0); return k+1
```

where `s` is at position `l` in `V<:MultiValue`.
"""
function _comp_wize_set_value!(r::AbstractMatrix{V},i,s::T,k,l) where {V,T}
  z = zero(T)
  ncomp = num_indep_components(V)
  r[i,k] = ntuple(i -> ifelse(i == l, s, z),Val(ncomp))
  return k + 1
end

function _gradient_nd!(
  b::CompWiseTensorPolyBasis{D,V,PT}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  s::MVector{D,T}, ::Val{K}) where {D,V,PT,G,T,K}

  if isHierarchical(PT)
    for d in 1:D
      _derivatives_1d!(PT,K,(c,g),x,d)
    end
  end

  k = 1
  @inbounds for (l,terms) in enumerate(get_comp_terms(b))

    if !isHierarchical(PT)
      for d in 1:D
        kld = b.orders[l,d]
        if isone(l) || kld ≠ b.orders[l-1,d]
          _derivatives_1d!(PT,kld,(c,g),x,d)
        end
      end
    end

    for ci in terms
      s .= one(T)

      for q in 1:D
        for d in 1:D
          if d != q
            s[q] *= c[d,ci[d]]
          else
            s[q] *= g[d,ci[d]]
          end
        end
      end

      k = _comp_wize_set_derivative!(r,i,s,k,Val(l),V)
    end
  end
end

"""
    _comp_wize_set_derivative!(r::AbstractMatrix{G},i,s,k,::Val{l},::Type{V})

```
z = zero(s)
r[i,k]     = G(z…, ..., z…, s…, z…, ..., z…) = (Dbᵏ)(xi)
return k+1
```

where `s…` is the `l`ᵗʰ set of components. This is the gradient or hessian of
the `k`ᵗʰ basis polynomial, whose nonzero component in `V` is the `l`ᵗʰ.
"""
@generated function _comp_wize_set_derivative!(
  r::AbstractMatrix{G},i,s,k,::Val{l},::Type{V}) where {G,l,V}

  N_val_dims = length(size(V))
  s_size = size(G)[1:end-N_val_dims]

  body = "T = eltype(s); z = zero(T);"
  m = Array{String}(undef, size(G))
  m .= "z"

  for ci in CartesianIndices(s_size)
    m[ci,l] = "(@inbounds s[$ci])"
  end
  body *= "@inbounds r[i,k] = ($(join(tuple(m...), ", ")));"

  body = Meta.parse(string("begin ",body," end"))
  return Expr(:block, body ,:(return k+1))
end

# See _cartprod_set_derivative!(r::AbstractMatrix{G},i,s,k,::Type{V}) where {G,V<:AbstractSymTensorValue{D}} where D
@generated function _comp_wize_set_derivative!(
  r::AbstractMatrix{G},i,s,k,::Type{V}) where {G,V<:AbstractSymTensorValue{D}} where D

  @notimplemented
end

function _hessian_nd!(
  b::CompWiseTensorPolyBasis{D,V,PT}, x,
  r::AbstractMatrix{H}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  h::AbstractMatrix{T},
  s::MMatrix{D,D,T}, ::Val{K}) where {D,V,PT,H,T,K}

  if isHierarchical(PT)
    for d in 1:D
      _derivatives_1d!(PT,K,(c,g,h),x,d)
    end
  end


  k = 1
  @inbounds for (l,terms) in enumerate(get_comp_terms(b))

    if !isHierarchical(PT)
      for d in 1:D
        kld = b.orders[l,d]
        if isone(l) || kld ≠ b.orders[l-1,d]
          _derivatives_1d!(PT,kld,(c,g),x,d)
        end
      end
    end

    for ci in terms
      s .= one(T)

      for r in 1:D
        for q in 1:D
          for d in 1:D
            if d != q && d != r
              s[r,q] *= c[d,ci[d]]
            elseif d == q && d ==r
              s[r,q] *= h[d,ci[d]]
            else
              s[r,q] *= g[d,ci[d]]
            end
          end
        end
      end

      k = _comp_wize_set_derivative!(r,i,s,k,Val(l),V)
    end
  end
end

