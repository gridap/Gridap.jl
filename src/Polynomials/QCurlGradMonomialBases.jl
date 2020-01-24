
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
"""
function QCurlGradMonomialBasis{D}(::Type{T},order::Int) where {D,T}
  @assert T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  _t = tfill(order,Val{D-1}())
  t = (order+1,_t...)
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
