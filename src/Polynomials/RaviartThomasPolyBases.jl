"""
    RaviartThomasPolyBasis{D,V,PT} <: PolynomialBasis{D,V,PT}

Basis of the vector valued (`V<:VectorValue{D}`) space

𝓡𝓣ᴰₙ = (𝓢ₙ)ᴰ ⊕ x (𝓢ₙ\\𝓢₍ₙ₋₁₎)

where 𝓢ₙ is a `D`-multivariate scalar polynomial space of maximum degree n = `K`-1.

This 𝓡𝓣ᴰₙ is the polynomial space for Raviart-Thomas elements with divergence in 𝓢ₙ.
Its maximum degree, that `get_order` returns, is n+1 = `K`.

!!! warning
    Using this basis on simplices is not recommanded, [`BarycentricPmΛBasis`](@ref) is better numerically conditioned for higher degrees, they are obtained by using `Bernstein` as argument of [`FEEC_poly_basis`](@ref) .

# Example:

```@example
# a basis for Raviart-Thomas on tetrahedra with divergence in 𝓟₂
b = RaviartThomasPolyBasis{3}(Monomial, Float64, 2)

# a basis for Raviart-Thomas on quadrilateral with divergence in 𝓟₃
b = RaviartThomasPolyBasis{2}(Monomial, Float64, 3, _q_filter)
```

The space 𝓢ₙ, typically 𝓟ᴰₙ or 𝓠ᴰₙ, does not need to have a tensor product
structure of 1D scalar spaces. Thus, the 𝓡𝓣ᴰₙ component's scalar spaces are not
tensor products either.

𝓢ₙ is defined like a scalar valued [`CartProdPolyBasis`](@ref) via the `_filter`
argument of the constructor, by default `_p_filter` for 𝓟ᴰₙ.
As a consequence, `PT` must be hierarchical, see [`isHierarchical`](@ref).
"""
struct RaviartThomasPolyBasis{D,V,PT} <: PolynomialBasis{D,V,PT}
  max_order::Int
  pterms::Vector{CartesianIndex{D}}
  sterms::Vector{CartesianIndex{D}}

  @doc"""
      RaviartThomasPolyBasis{D}(::Type{PT}, ::Type{T}, order::Int, _filter::Function=_p_filter)

  Where `_filter` defines 𝓢ₙ and `order` = n = K-1 (cf. struct docstring).
  """
  function RaviartThomasPolyBasis{D}(
    ::Type{PT}, ::Type{T}, order::Int,
    _filter::Function=_p_filter
    ) where {PT<:Polynomial,D,T}

    @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
    @check D > 1
    @check isconcretetype(PT) "PT needs to be a concrete <:Polynomial type"
    @check isHierarchical(PT) "The polynomial basis must be hierarchical for this space."

    V = VectorValue{D,T}
    indexbase = 1

    # terms defining 𝓢ₙ
    P_k = MonomialBasis(Val(D), T, order, _filter)
    pterms = P_k.terms
    msg =  "Some term defining `𝓢ₙ` contain a higher index than the maximum,
      `order`+1, please fix the `_filter` argument"
    @check all( pterm -> (maximum(Tuple(pterm) .- indexbase, init=0) <= order), pterms) msg

    # terms defining 𝓢ₙ\𝓢ₙ₋₁
    _minus_one_order_filter = term -> _filter(Tuple(term) .- indexbase, order-1)
    sterms = filter(!_minus_one_order_filter, pterms)

    new{D,V,PT}(order+1,pterms,sterms)
  end
end

Base.size(b::RaviartThomasPolyBasis{D}) where {D} = (D*length(b.pterms) + length(b.sterms), )
get_order(b::RaviartThomasPolyBasis) = b.max_order
get_orders(b::RaviartThomasPolyBasis{D}) where D = tfill(get_order(b), Val(D))

function testvalue(::Type{RaviartThomasPolyBasis{D,V,PT}}) where {D,V,PT}
  T = eltype(V)
  RaviartThomasPolyBasis{D}(PT, T, 0)
end

#################################
# nD evaluations implementation #
#################################

function _evaluate_nd!(
  b::RaviartThomasPolyBasis{D,V,PT}, x,
  r::AbstractMatrix{Vr}, i,
  c::AbstractMatrix{T}, K) where {D,V,PT,Vr,T}

  for d in 1:D
    _evaluate_1d!(PT,K,c,x,d)
  end

  m = zero(Mutable(Vr))
  k = 1

  @inbounds begin
    for l in 1:D
      for ci in b.pterms

        s = one(T)
        for d in 1:D
          s *= c[d,ci[d]]
        end

        k = _comp_wize_set_value!(r,i,s,k,l)
      end
    end

    for ci in b.sterms
      for i in 1:D
        m[i] = zero(T)
      end

      for l in 1:D

        s = x[l]
        for d in 1:D
          s *= c[d,ci[d]]
        end

        m[l] = s
      end

      r[i,k] = m
      k += 1
    end

  end
end

function _gradient_nd!(
  b::RaviartThomasPolyBasis{D,V,PT}, x,
  r::AbstractMatrix{G}, i,
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  s::MVector{D,T}, K) where {D,V,PT,G,T}

  for d in 1:D
    _derivatives_1d!(PT,K,(c,g),x,d)
  end

  m = zero(Mutable(G))
  k = 1

  @inbounds begin
    for l in 1:D
      for ci in b.pterms

        for i in eachindex(s)
          s[i] = one(T)
        end

        for q in 1:D
          for d in 1:D
            if d != q
               s[q] *= c[d,ci[d]]
            else
               s[q] *= g[d,ci[d]]
            end
          end
        end

        k = _comp_wize_set_derivative!(r,i,s,k,l,V)
      end
    end

    for ci in b.sterms

      for i in eachindex(m)
        m[i] = zero(T)
      end

      for l in 1:D

        for i in eachindex(s)
          s[i] = x[l]
        end

        aux = one(T)
        for q in 1:D
          aux *= c[q,ci[q]]
          for d in 1:D
            if d != q
               s[q] *= c[d,ci[d]]
            else
               s[q] *= g[d,ci[d]]
            end
          end
        end
        s[l] += aux

        for i in 1:D
           m[i,l] = s[i]
        end
      end
      r[i,k] = m
      k += 1
    end

  end
end

