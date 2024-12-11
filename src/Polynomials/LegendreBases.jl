"""
    Legendre <: Polynomial

Type representing the Legendre polynomials
"""
struct Legendre <: Polynomial   end

isHierarchical(::Legendre) = true

"""
    LegendreBasis{D,T,K} = TensorPolynomialBasis{D,T,K,Legendre}

Multivariate scalar' or `Multivalue`'d Legendre basis, see [`TensorPolynomialBasis`](@ref)
"""
const LegendreBasis{D,T,K} = TensorPolynomialBasis{D,T,K,Legendre}

"""
    LegendreBasis{D}(::Type{T}, order::Int, terms::Vector) where {D,T}
    LegendreBasis{D}(::Type{T}, orders::Tuple [, filter::Function]) where {D,T}
    LegendreBasis{D}(::Type{T}, order::Int [, filter::Function]) where {D,T}

Convenience constructors of LegendreBasis{D,T}.
"""
LegendreBasis{D}(args...) where {D} = TensorPolynomialBasis{D}(Legendre, args...)

# 1D evaluation implementation

function _evaluate_1d!(::Type{Legendre}, k, v::AbstractMatrix{T},x,d) where T<:Number
  n = k + 1
  o = one(T)
  @inbounds v[d,1] = o
  if n > 1
    ξ = ( 2*x[d] - 1 )
    for i in 2:n
      # The sqrt(2i-1) factor normalizes the basis polynomial for L2 scalar product on ξ∈[0,1], indeed:
      # ∫[0,1] Pn(2ξ-1)^2 dξ = 1/2 ∫[-1,1] Pn(t)^2 dt = 1/(2n+1)
      # C.f. Eq. (1.25) in Section 1.1.5 in Ern & Guermond book (2013).
      @inbounds v[d,i] = sqrt(2*i-1)*jacobi(ξ,i-1,0,0)
    end
  end
end

function _gradient_1d!(::Type{Legendre}, k, g::AbstractMatrix{T},x,d) where T<:Number
  n = k + 1
  z = zero(T)
  @inbounds g[d,1] = z
  if n > 1
    ξ = ( 2*x[d] - 1 )
    for i in 2:n
      @inbounds g[d,i] = sqrt(2*i-1)*i*jacobi(ξ,i-2,1,1)
    end
  end
end

function _hessian_1d!(::Type{Legendre}, k, h::AbstractMatrix{T},x,d) where T<:Number
  n = k + 1
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
