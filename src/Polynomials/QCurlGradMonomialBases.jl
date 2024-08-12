
"""
    struct QCurlGradMonomialBasis{...} <: AbstractArray{Monomial}

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Raviart-Thomas reference elements on n-cubes.
The type parameters and fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct QCurlGradMonomialBasis{D,T} <: AbstractVector{Monomial}
  qgrad::QGradMonomialBasis{D,T}
  function QCurlGradMonomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    qgrad = QGradMonomialBasis(T,order,terms,perms)
    new{D,T}(qgrad)
  end
end

Base.size(a::QCurlGradMonomialBasis) = (length(a.qgrad),)
# @santiagobadia : Not sure we want to create the monomial machinery
Base.getindex(a::QCurlGradMonomialBasis,i::Integer) = Monomial()
Base.IndexStyle(::QCurlGradMonomialBasis) = IndexLinear()

"""
    QCurlGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}

Returns a `QCurlGradMonomialBasis` object. `D` is the dimension
of the coordinate space and `T` is the type of the components in the vector-value.
The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the Q space of degree `order`.
"""
function QCurlGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  _order = order+1
  _t = tfill(_order,Val{D-1}())
  t = (_order+1,_t...)
  terms = CartesianIndices(t)
  perms = _prepare_perms(D)
  QCurlGradMonomialBasis(T,order,terms,perms)
end

# @santiagobadia: This is dirty, I would put here VectorValue{D,T}
return_type(::QCurlGradMonomialBasis{D,T}) where {D,T} = T

function return_cache(f::QCurlGradMonomialBasis,x::AbstractVector{<:Point})
  return_cache(f.qgrad,x)
end

function evaluate!(cache,f::QCurlGradMonomialBasis,x::AbstractVector{<:Point})
  evaluate!(cache,f.qgrad,x)
end

function return_cache(
  fg::FieldGradientArray{N,<:QCurlGradMonomialBasis},
  x::AbstractVector{<:Point}) where N

  f = fg.fa
  return_cache(FieldGradientArray{N}(f.qgrad),x)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{N,<:QCurlGradMonomialBasis},
  x::AbstractVector{<:Point}) where N

  f = fg.fa
  evaluate!(cache,FieldGradientArray{N}(f.qgrad),x)
end

"""
    num_terms(f::QCurlGradMonomialBasis{D,T}) where {D,T}
"""
num_terms(f::QCurlGradMonomialBasis{D,T}) where {D,T} = length(f.qgrad.terms)*D

get_order(f::QCurlGradMonomialBasis{D,T}) where {D,T} = get_order(f.qgrad)
