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
# @santiagobadia : I have added field to have control whether k returns a field
# or not. If it does not return a field, just call composition
# We could create FieldOperation an eliminate field, probably too much...

@inline field_composition(k,l...) = GenericField(composition(k,l...))
@inline field_array_composition(k,l...) = GenericFieldArray(composition(k,l...))

# @inline function field_composition(k,l::Vararg{<:NewField})
  # FieldComposition(k,l...)
# end

# """
#     apply_mapping_gradient(k,f...)

# Returns a field representing the gradient of the field obtained with

#     apply_mapping_to_field(k,f...)
# """
# struct FieldComposition{K,L} <: NewField
#   k::K
#   l::L
#   @inline function FieldComposition(k,l...)
#     new{typeof(k),typeof(l)}(k,l)
#   end
# end

# function return_type(c::FieldComposition,x)
#   Ts = return_types(c.l,x)
#   return_type(c.k, testvalues(Ts...)...)
# end

# function return_cache(c::FieldComposition,x)
#   cl = return_caches(c.l,x)
#   lx = evaluate!(cl,c.l,x)
#   ck = return_cache(c.k,lx...)
#   (ck,cl)
# end

# @inline function evaluate!(cache,c::FieldComposition,x)
#   ck, cf = cache
#   lx = evaluate!(cf,c.l,x)
#   evaluate!(ck,c.k,lx...)
# end

# @inline function field_composition(k,l::Vararg{<:NewField})
#   FieldComposition(k,l...)
# end

# """
#     apply_mapping_gradient(k,f...)

# Returns a field representing the gradient of the field obtained with

#     apply_mapping_to_field(k,f...)
# """
# struct FieldComposition{K,L} <: NewField
#   k::K
#   l::L
#   @inline function FieldComposition(k,l...)
#     new{typeof(k),typeof(l)}(k,l)
#   end
# end

# function return_type(c::FieldComposition,x)
#   Ts = return_types(c.l,x)
#   return_type(c.k, testvalues(Ts...)...)
# end

# function return_cache(c::FieldComposition,x)
#   cl = return_caches(c.l,x)
#   lx = evaluate!(cl,c.l,x)
#   ck = return_cache(c.k,lx...)
#   (ck,cl)
# end

# @inline function evaluate!(cache,c::FieldComposition,x)
#   ck, cf = cache
#   lx = evaluate!(cf,c.l,x)
#   evaluate!(ck,c.k,lx...)
# end
