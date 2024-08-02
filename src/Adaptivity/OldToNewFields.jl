
abstract type NewFieldType end;
struct CoarsenedNewFieldType <: NewFieldType end;
struct RefinedOrUntouchedNewFieldType <: NewFieldType end;   

# Unfortunately, I cannot use traits in OldToNewField as its 
# usage results in conversion errors when leveraging it with 
# the lazy_map infrastructure
struct OldToNewField <: Fields.Field
  new_field_type::NewFieldType
  fine_to_coarse_field
  refined_or_untouched_field
  function OldToNewField(a::NewFieldType,b::Fields.Field,c::Fields.Field)
    new(a,b,c)
  end
end

function OldToNewField(old_fields::AbstractArray{<:Fields.Field},
                       rrule::RefinementRule,
                       child_ids::AbstractVector{<:Integer})
  @check length(old_fields)==length(child_ids)                   
  if length(old_fields) == 1
    cell_map  = get_cell_map(rrule)[child_ids[1]]
    old_field = old_fields[1]
    fine_to_coarse_field = FineToCoarseField(
          [old_field for i = 1:num_subcells(rrule)],
          rrule,
          [i for i=1:num_subcells(rrule)])
    refined_or_untouched_field = old_field∘cell_map
    return OldToNewField(RefinedOrUntouchedNewFieldType(),fine_to_coarse_field,refined_or_untouched_field)
  else 
    @check length(old_fields) <= num_subcells(rrule)
    fine_to_coarse_field = FineToCoarseField(old_fields,rrule,child_ids)
    cell_map = get_cell_map(rrule)[1]
    refined_or_untouched_field = old_fields[1]∘cell_map 
    return OldToNewField(CoarsenedNewFieldType(),fine_to_coarse_field,refined_or_untouched_field)
  end
end

function Fields.return_cache(a::OldToNewField,x::AbstractArray{<:Point})
  f2c_cache = Fields.return_cache(a.fine_to_coarse_field,x)
  rou_cache = Fields.return_cache(a.refined_or_untouched_field,x)
  return (f2c_cache,rou_cache)
end

function Fields.evaluate!(cache,a::OldToNewField,x::AbstractArray{<:Point})
  if isa(a.new_field_type,CoarsenedNewFieldType)
    f2c_cache,rou_cache = cache
    Fields.evaluate!(f2c_cache,a.fine_to_coarse_field,x)
  else
    @assert isa(a.new_field_type,RefinedOrUntouchedNewFieldType)
    f2c_cache,rou_cache = cache
    Fields.evaluate!(rou_cache,a.refined_or_untouched_field,x)
  end  
end
