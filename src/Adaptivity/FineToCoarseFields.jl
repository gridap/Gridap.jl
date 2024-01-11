
"""
  Given a domain and a non-overlapping refined cover, a `FineToCoarseField`
  is a `Field` defined in the domain and constructed by a set of fields defined on 
  the subparts of the covering partition.
  The refined cover is represented by a `RefinementRule`. 
"""
struct FineToCoarseField{A<:AbstractArray{<:Field},B<:RefinementRule} <: Field
  fine_fields :: A
  rrule       :: B
  is_zero     :: Vector{Bool}
  function FineToCoarseField(
    fine_fields::AbstractArray{<:Field},
    rrule::RefinementRule,
    is_zero::Vector{Bool}=fill(false,num_subcells(rrule))
  )
    @check length(fine_fields) == num_subcells(rrule)
    A = typeof(fine_fields)
    B = typeof(rrule)
    new{A,B}(fine_fields,rrule,is_zero)
  end
end

# Necessary for distributed meshes, where not all children of a coarse cell may belong to the processor. 
function FineToCoarseField(
  fine_fields::AbstractArray{T},
  rrule::RefinementRule,
  child_ids::AbstractArray{<:Integer}
) where T <: Field
  fields  = Vector{Union{T,ZeroField{T}}}(undef,num_subcells(rrule))
  fields  = fill!(fields,ZeroField(testitem(fine_fields)))
  is_zero = fill(true,num_subcells(rrule))
  for (k,id) in enumerate(child_ids)
    fields[id] = fine_fields[k]
    is_zero[id] = false
  end
  return FineToCoarseField(fields,rrule,is_zero)
end

function Geometry.return_cache(a::FineToCoarseField,x::AbstractArray{<:Point})
  fields, rr, is_zero = a.fine_fields, a.rrule, a.is_zero
  cmaps = get_inverse_cell_map(a.rrule)

  xi_cache = array_cache(x)
  fi_cache = array_cache(fields)
  mi_cache = array_cache(cmaps)

  # Generic caches
  child_ids = map(i -> x_to_cell(rr,getindex!(xi_cache,x,i)),eachindex(x))
  pos = findfirst(id -> !is_zero[id],child_ids)
  xi = getindex!(xi_cache,x,pos)
  id = child_ids[pos]
  id = x_to_cell(rr,xi)
  mi = getindex!(mi_cache,cmaps,id)
  fi = getindex!(fi_cache,fields,id)

  zi_cache = Fields.return_cache(mi,xi)
  zi = evaluate!(zi_cache,mi,xi)

  yi_type  = Fields.return_type(fi,zi)
  y_cache  = Arrays.CachedArray(zeros(yi_type,size(x)))

  # Evaluation caches
  fi_zero = ZeroField(fi)
  yi_nonzero_cache = Fields.return_cache(fi,zi)
  yi_zero_cache = Fields.return_cache(fi_zero,zi)
  yi_cache = (yi_nonzero_cache,yi_zero_cache)

  return fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache
end

function Geometry.evaluate!(cache,a::FineToCoarseField,x::AbstractArray{<:Point})
  fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache = cache
  fields, rr, is_zero = a.fine_fields, a.rrule, a.is_zero
  cmaps = get_inverse_cell_map(rr)

  Arrays.setsize!(y_cache, size(x))

  for i in eachindex(x)
    xi = getindex!(xi_cache,x,i)
    child_id = x_to_cell(rr,xi)
    fi = getindex!(fi_cache,fields,child_id)
    mi = getindex!(mi_cache,cmaps,child_id)
    zi = Fields.evaluate!(zi_cache,mi,xi)
    _yi_cache = yi_cache[is_zero[child_id]+1]
    y_cache.array[i] = Fields.evaluate!(_yi_cache,fi,zi)
  end
  return y_cache.array
end

# Fast evaluation of FineToCoarseFields: 
# Points are pre-classified into the children cells, which allows for the search to be 
# skipped entirely. 
function Geometry.return_cache(a::FineToCoarseField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  fields, is_zero = a.fine_fields, a.is_zero
  cmaps = get_inverse_cell_map(a.rrule)

  # Generic caches
  xi_cache = array_cache(x)
  fi_cache = array_cache(fields)
  mi_cache = array_cache(cmaps)
  id_cache = array_cache(child_ids)

  pos = findfirst(id -> !is_zero[id],child_ids)
  xi = getindex!(xi_cache,x,pos)
  id = getindex!(id_cache,child_ids,pos)
  mi = getindex!(mi_cache,cmaps,id)
  fi = getindex!(fi_cache,fields,id)

  zi_cache = Fields.return_cache(mi,xi)
  zi = evaluate!(zi_cache,mi,xi)

  yi_type  = Fields.return_type(fi,zi)
  y_cache  = Arrays.CachedArray(zeros(yi_type,size(x)))

  # Evaluation caches
  fi_zero = ZeroField(fi)
  yi_nonzero_cache = Fields.return_cache(fi,zi)
  yi_zero_cache = Fields.return_cache(fi_zero,zi)
  yi_cache = (yi_nonzero_cache,yi_zero_cache)

  return fi_cache, mi_cache, xi_cache, id_cache, zi_cache, yi_cache, y_cache
end

function Geometry.evaluate!(cache,a::FineToCoarseField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  fi_cache, mi_cache, xi_cache, id_cache, zi_cache, yi_cache, y_cache = cache
  cmaps  = get_inverse_cell_map(a.rrule)
  fields, is_zero = a.fine_fields, a.is_zero

  Arrays.setsize!(y_cache, size(x))

  for i in eachindex(x)
    xi = getindex!(xi_cache,x,i)
    id = getindex!(id_cache,child_ids,i)
    fi = getindex!(fi_cache,fields,id)
    mi = getindex!(mi_cache,cmaps,id)
    zi = Fields.evaluate!(zi_cache,mi,xi)
    _yi_cache = yi_cache[is_zero[id]+1]
    y_cache.array[i] = Fields.evaluate!(_yi_cache,fi,zi)
  end
  return y_cache.array
end
