"""
    Chebyshev{kind} <: Polynomial

Type representing Chebyshev polynomials of the
- first kind: `Chebyshev{:T}`
- second kind: `Chebyshev{:U}`
C.f. [Chebyshev polynomials](@ref) section.
"""
struct Chebyshev{kind} <: Polynomial
  Chebyshev{:T}() = new{:T}()
  Chebyshev{:U}() = new{:U}()
end

isHierarchical(::Type{<:Chebyshev}) = true

"""
    ChebyshevBasis{D,V,kind} = CartProdPolyBasis{D,V,Chebyshev{kind}}

Alias for cartesian product Chebyshev basis, scalar valued or multivalued.
"""
const ChebyshevBasis{D,V,kind} = CartProdPolyBasis{D,V,Chebyshev{kind}}

"""
    ChebyshevBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector; kind=:T)
    ChebyshevBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function; kind=:T])
    ChebyshevBasis(::Val{D}, ::Type{V}, orders::Tuple [, filter::Function; kind=:T])

High level constructors of [`ChebyshevBasis`](@ref).
"""
ChebyshevBasis(args...; kind=:T) = CartProdPolyBasis(Chebyshev{kind}, args...)

function CartProdPolyBasis(
  ::Type{Chebyshev{:U}}, ::Val{D}, ::Type{V}, ::Int) where {D, V}

  @notimplemented "1D evaluation for second kind need to be implemented here"
end


# 1D evaluation implementation

function _evaluate_1d!(
  ::Type{Chebyshev{kind}},K,c::AbstractMatrix{T},x,d) where {kind,T<:Number}

  if iszero(K)
    @inbounds c[1] = one(T)
    return
  end

  n = K + 1        # n > 1
  ξ = (2*x[d] - 1) # ξ ∈ [-1,1]
  ξ2 = 2*ξ

  @inbounds c[d,1] = one(T)
  @inbounds c[d,2] = (kind == :T) ? ξ :  ξ2
  for i in 3:n
    @inbounds c[d,i] = c[d,i-1]*ξ2 - c[d,i-2]
  end
end

function _gradient_1d!(
  ::Type{Chebyshev{:T}},K,g::AbstractMatrix{T},x,d) where T<:Number

  if iszero(K)
    @inbounds g[d,1] = zero(T)
    return
  end

  n = K + 1  # n>1
  z = zero(T)
  o = one(T)
  ξ = T(2*x[d] - 1)
  dξdx = T(2.0)

  unm1 = o
  un = 2*ξ
  @inbounds g[d,1] = z               # dT_0 = 0
  @inbounds g[d,2] = dξdx*o          # dT_1 = 1*U_0 = 1
  for i in 3:n
    @inbounds g[d,i] = dξdx*(i-1)*un # dT_i = i*U_{i-1}
    un, unm1 = 2*ξ*un - unm1, un
  end
end

_evaluate_1d!(::Type{Chebyshev{:U}},K,h::AbstractMatrix{T},x,d) where T<:Number = @notimplemented
_gradient_1d!(::Type{Chebyshev{:U}},K,h::AbstractMatrix{T},x,d) where T<:Number = @notimplemented
_hessian_1d!( ::Type{<:Chebyshev},  K,h::AbstractMatrix{T},x,d) where T<:Number = @notimplemented

