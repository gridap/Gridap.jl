
"""
    struct QCurlGradMonomialBasis{...} <: Field

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Raviart-Thomas reference elements on n-cubes.
The type parameters and fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct QCurlGradMonomialBasis{D,T} <: Field
  qgrad::QGradMonomialBasis{D,T}
  function QCurlGradMonomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    qgrad = QGradMonomialBasis(T,order,terms,perms)
    new{D,T}(qgrad)
  end
end

"""
    QCurlGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}

Returns a `QCurlGradMonomialBasis` object. `D` is the dimension
of the coordinate space and `T` is the type of the components in the vector-value.
The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the Q space of degree `order`.
"""
function QCurlGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}
  @assert T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  _order = order+1
  _t = tfill(_order,Val{D-1}())
  t = (_order+1,_t...)
  terms = CartesianIndices(t)
  perms = _prepare_perms(D)
  QCurlGradMonomialBasis(T,order,terms,perms)
end

function field_cache(f::QCurlGradMonomialBasis,x)
  field_cache(f.qgrad,x)
end

@inline function evaluate_field!(cache,f::QCurlGradMonomialBasis,x)
  evaluate_field!(cache,f.qgrad,x)
end

function gradient_cache(f::QCurlGradMonomialBasis,x)
  gradient_cache(f.qgrad,x)
end

@inline function evaluate_gradient!(cache,f::QCurlGradMonomialBasis,x)
  evaluate_gradient!(cache,f.qgrad,x)
end

get_value_type(::QCurlGradMonomialBasis{D,T}) where {D,T} = T

"""
    num_terms(f::QCurlGradMonomialBasis{D,T}) where {D,T}
"""
num_terms(f::QCurlGradMonomialBasis{D,T}) where {D,T} = length(f.qgrad.terms)*D

get_order(f::QCurlGradMonomialBasis{D,T}) where {D,T} = get_order(f.qgrad)

