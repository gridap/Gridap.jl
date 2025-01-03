"""
    Chebyshev{kind} <: Polynomial

Type representing Chebyshev polynomials of first and second kind where `kind` is
either `:T` or `:U` for first and second kind Chebyshev polynomials respectively.
"""
struct Chebyshev{kind} <: Polynomial end

isHierarchical(::Type{<:Chebyshev}) = true

"""
    ChebyshevBasis{D,V,kind,K} = UniformPolyBasis{D,V,K,Chebyshev{kind}}

Alias for Chebyshev multivariate scalar' or `Multivalue`'d basis, see [`UniformPolyBasis`](@ref).
"""
const ChebyshevBasis{D,V,kind,K} = UniformPolyBasis{D,V,K,Chebyshev{kind}}

"""
    ChebyshevBasis(::Val{D}, ::Type{V}, order::Int, terms::Vector; kind=:T)
    ChebyshevBasis(::Val{D}, ::Type{V}, order::Int [, filter::Function; kind=:T])
    ChebyshevBasis(::Val{D}, ::Type{V}, orders::Tuple [, filter::Function; kind=:T])

High level constructors of [`ChebyshevBasis`](@ref).
"""
ChebyshevBasis(args...; kind=:T) = UniformPolyBasis(Chebyshev{kind}, args...)

UniformPolyBasis{D}(::Type{Chebyshev{:U}}, args...) where D = @notimplemented "1D evaluation for second kind needed here"


# 1D evaluation implementation

function _evaluate_1d!(
  ::Type{<:Chebyshev},::Val{0},c::AbstractMatrix{T},x,d) where T<:Number

  @inbounds c[d,1] = one(T)
end

function _evaluate_1d!(
  ::Type{Chebyshev{:T}},::Val{K},c::AbstractMatrix{T},x,d) where {K,T<:Number}

  n = K + 1        # n > 1
  ξ = (2*x[d] - 1) # ξ ∈ [-1,1]
  ξ2 = 2*ξ

  @inbounds c[d,1] = one(T)
  @inbounds c[d,2] = ξ
  for i in 3:n
    @inbounds c[d,i] = c[d,i-1]*ξ2 - c[d,i-2]
  end
end


function _gradient_1d!(
  ::Type{<:Chebyshev},::Val{0},g::AbstractMatrix{T},x,d) where T<:Number

  @inbounds g[d,1] = zero(T)
end

function _gradient_1d!(
  ::Type{Chebyshev{:T}},::Val{K},g::AbstractMatrix{T},x,d) where {K,T<:Number}

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


function _hessian_1d!(
  ::Type{Chebyshev{:T}},::Val{K},h::AbstractMatrix{T},x,d) where {K,T<:Number}

  @notimplemented
end


_evaluate_1d!(::Type{Chebyshev{:U}},::Val{K},c::AbstractMatrix{T},x,d) where {K,T<:Number} = @notimplemented
_gradient_1d!(::Type{Chebyshev{:U}},::Val{K},g::AbstractMatrix{T},x,d) where {K,T<:Number} = @notimplemented
_hessian_1d!( ::Type{Chebyshev{:U}},::Val{K},h::AbstractMatrix{T},x,d) where {K,T<:Number} = @notimplemented
