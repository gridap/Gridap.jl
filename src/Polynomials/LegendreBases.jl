"""
    Legendre <: Polynomial

Type representing the Legendre polynomials.
"""
struct Legendre <: Polynomial   end

isHierarchical(::Type{Legendre}) = true

"""
    LegendreBasis{D,V,K} = UniformPolyBasis{D,V,K,Legendre}

Alias for Legendre multivariate scalar' or `Multivalue`'d basis, see [`UniformPolyBasis`](@ref).
"""
const LegendreBasis{D,V,K} = UniformPolyBasis{D,V,K,Legendre}

"""
    LegendreBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector)
    LegendreBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function])
    LegendreBasis(::Val{D}, ::Type{V}, orders::Tuple [, filter::Function])

High level constructors of [`LegendreBasis`](@ref).
"""
LegendreBasis(args...) = UniformPolyBasis(Legendre, args...)


# 1D evaluation implementation

# TODO optimize evaluation by using the iterative formula explicitely

function _evaluate_1d!(::Type{Legendre},::Val{K},c::AbstractMatrix{T},x,d) where {K,T<:Number}
  n = K + 1
  @inbounds c[d,1] = one(T)
  if n > 1
    ξ = ( 2*x[d] - 1 )
    for i in 2:n
      # The sqrt(2i-1) factor normalizes the basis polynomial for L2 scalar
      # product on ξ∈[0,1], indeed:
      # ∫[0,1] Pn(2ξ-1)^2 dξ = 1/2 ∫[-1,1] Pn(t)^2 dt = 1/(2n+1)
      # C.f. Eq. (1.25) in Section 1.1.5 in Ern & Guermond book (2013).
      @inbounds c[d,i] = sqrt(2*i-1)*jacobi(ξ,i-1,0,0)
    end
  end
end

function _gradient_1d!(::Type{Legendre},::Val{K},g::AbstractMatrix{T},x,d) where {K,T<:Number}
  n = K + 1
  z = zero(T)
  @inbounds g[d,1] = z
  if n > 1
    ξ = ( 2*x[d] - 1 )
    for i in 2:n
      @inbounds g[d,i] = sqrt(2*i-1)*i*jacobi(ξ,i-2,1,1)
    end
  end
end

function _hessian_1d!(::Type{Legendre},::Val{K},h::AbstractMatrix{T},x,d) where {K,T<:Number}
  n = K + 1
  z = zero(T)
  @inbounds h[d,1] = z
  if n > 1
    @inbounds h[d,2] = z
    ξ = ( 2*x[d] - 1 )
    for i in 3:n
      @inbounds h[d,i] = sqrt(2*i-1)*(i*(i+1)/2)*jacobi(ξ,i-3,2,2)
    end
  end
end
