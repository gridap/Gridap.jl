"""
    apply_mapping_to_field(k,f...) -> Field

Returns a field obtained by applying the mapping `k` to the
values of the fields in `f`. That is, the returned field evaluated at
a vector of points `x` provides the value obtained by applying mapping `k` to the
values of the fields `f` at the vector of points `x`. Formally, the resulting field at a
vector of points
 `x` is defined as

    fx = evaluate_fields(f,x)
    apply_mapping(k,fx...)

In order to be able to call the [`field_gradient`](@ref) function of the
resulting field, one needs to define the gradient operator
associated with the underlying mapping.
This is done by adding a new method to [`apply_mapping_gradient(k,f...)`](@ref) for each mapping type.
"""
# @inline function composition(k::Mapping,l::NTuple{N,<:Mapping}) where N
# @santiagobadia : To decide tuple or no tuple
# @inline function composition(k::Mapping,l::Tuple{Vararg{<:Mapping}})
@inline function composition(k,l...)
  CompositionMapping(k,l...)
end

# @inline function composition(k,l)
#   CompositionMapping(k,l)
# end

"""
    apply_mapping_gradient(k,f...)

Returns a field representing the gradient of the field obtained with

    apply_mapping_to_field(k,f...)
"""
# function apply_mapping_gradient(k,f...)
#   @abstractmethod
# end

# @inline apply_mapping_gradient(k::BCasted{typeof(+)},a) = field_gradient(a)

# @inline apply_mapping_gradient(k::BCasted{typeof(-)},a...) = apply_mapping_to_field(k,field_gradients(a...)...)

# @inline apply_mapping_gradient(k::BCasted{typeof(+)},a...) = apply_mapping_to_field(k,field_gradients(a...)...)

# Result of applying a mapping to the value of some fields

struct CompositionMapping{K,L} <: Mapping
  k::K
  l::L
  @inline function CompositionMapping(k,l...)
    new{typeof(k),typeof(l)}(k,l)
  end
end

function return_type(c::CompositionMapping,x)
  Ts = return_types(c.l,x)
  return_type(c.k, testvalues(Ts...)...)
end

function return_cache(c::CompositionMapping,x)
  cl = return_caches(c.l,x)
  lx = evaluate!(cl,c.l,x)
  ck = return_cache(c.k,lx...)
  (ck,cl)
end

@inline function evaluate!(cache,c::CompositionMapping,x)
  ck, cf = cache
  lx = evaluate!(cf,c.l,x)
  evaluate!(ck,c.k,lx...)
end

# @inline function field_gradient(f::CompositionMapping)
  # apply_mapping_gradient(f.k,f.f...)
# end
