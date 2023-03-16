
"""
Given a domain and a non-overlapping refined cover, a `FineToCoarseField`
is a `Field` defined in the domain and constructed by a set of fields defined on 
the subparts of the covering partition.
The refined cover is represented by a `RefinementRule`. 
"""
struct FineToCoarseField{A<:AbstractArray{<:Field},B<:RefinementRule} <: Field
  fine_fields :: A
  rrule       :: B
  function FineToCoarseField(fine_fields::AbstractArray{<:Field},rrule::RefinementRule)
    @check length(fine_fields) == num_subcells(rrule)
    A = typeof(fine_fields)
    B = typeof(rrule)
    new{A,B}(fine_fields,rrule)
  end
end

# Necessary for distributed meshes, where not all children of a coarse cell may belong to the processor. 
function FineToCoarseField(fine_fields::AbstractArray{<:Field},rrule::RefinementRule,child_ids::AbstractArray{<:Integer})
  fields = Vector{Field}(undef,num_subcells(rrule))
  fields = fill!(fields,ConstantField(0.0))
  for (k,id) in enumerate(child_ids)
    fields[id] = fine_fields[k]
  end
  return FineToCoarseField(fields,rrule)
end

"""
function Geometry.return_cache(a::FineToCoarseField,x::Point)
  fields = a.fine_fields
  cmaps = get_inverse_cell_map(a.rrule)

  fi_cache = array_cache(fields)
  cm_cache = array_cache(cmaps)
  xi_cache = Fields.return_cache(first(cmaps),x)
  yi_cache = Fields.return_cache(first(fields),x)
  return fi_cache, cm_cache, xi_cache, yi_cache
end

function Geometry.evaluate!(cache,a::FineToCoarseField,x::Point)
  fi_cache, cm_cache, xi_cache, yi_cache = cache
  fields, rr = a.fine_fields, a.rrule
  cmaps = get_inverse_cell_map(rr)

  child_id = x_to_cell(rr,x) # Find correct subcell
  fi = getindex!(fi_cache,fields,child_id)
  mi = getindex!(cm_cache,cmaps,child_id)
  xi = Fields.evaluate!(xi_cache,mi,x)  # xc -> xf
  yi = Fields.evaluate!(yi_cache,fi,xi) # xf -> yf
  return yi
end
"""

function Geometry.return_cache(a::FineToCoarseField,x::AbstractArray{<:Point})
  fields, rr = a.fine_fields, a.rrule
  cmaps = get_inverse_cell_map(a.rrule)

  xi_cache = array_cache(x)
  fi_cache = array_cache(fields)
  mi_cache = array_cache(cmaps)

  xi = getindex!(xi_cache,x,1)
  child_id = x_to_cell(rr,xi)
  mi = getindex!(mi_cache,cmaps,child_id)
  fi = getindex!(fi_cache,fields,child_id)

  zi_cache = Fields.return_cache(mi,xi)
  zi = evaluate!(zi_cache,mi,xi)

  yi_type  = Fields.return_type(fi,zi)
  yi_cache = Fields.return_cache(fi,zi)
  y_cache  = Arrays.CachedArray(yi_type,1)

  return fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache
end

function Geometry.evaluate!(cache,a::FineToCoarseField,x::AbstractArray{<:Point})
  fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache = cache
  fields, rr = a.fine_fields, a.rrule
  cmaps = get_inverse_cell_map(rr)

  Arrays.setsize!(y_cache, size(x))

  for i in eachindex(x)
    xi = getindex!(xi_cache,x,i)
    child_id = x_to_cell(rr,xi)
    fi = getindex!(fi_cache,fields,child_id)
    mi = getindex!(mi_cache,cmaps,child_id)
    zi = Fields.evaluate!(zi_cache,mi,xi)
    y_cache.array[i] = Fields.evaluate!(yi_cache,fi,zi)
  end
  return y_cache.array
end

# Fast evaluation of FineToCoarseFields: 
# Points are pre-classified into the children cells, which allows for the search to be 
# skipped entirely. 
function Geometry.return_cache(a::FineToCoarseField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  fields = a.fine_fields
  cmaps = get_inverse_cell_map(a.rrule)

  xi_cache = array_cache(x)
  fi_cache = array_cache(fields)
  mi_cache = array_cache(cmaps)
  id_cache = array_cache(child_ids)

  xi = getindex!(xi_cache,x,1)
  id = getindex!(id_cache,child_ids,1)
  mi = getindex!(mi_cache,cmaps,id)
  fi = getindex!(fi_cache,fields,id)

  zi_cache = Fields.return_cache(mi,xi)
  zi = evaluate!(zi_cache,mi,xi)

  yi_type  = Fields.return_type(fi,zi)
  yi_cache = Fields.return_cache(fi,zi)
  y_cache  = Arrays.CachedArray(yi_type,1)

  return fi_cache, mi_cache, xi_cache, id_cache, zi_cache, yi_cache, y_cache
end

function Geometry.evaluate!(cache,a::FineToCoarseField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  fi_cache, mi_cache, xi_cache, id_cache, zi_cache, yi_cache, y_cache = cache
  cmaps  = get_inverse_cell_map(a.rrule)
  fields = a.fine_fields

  Arrays.setsize!(y_cache, size(x))

  for i in eachindex(x)
    xi = getindex!(xi_cache,x,i)
    id = getindex!(id_cache,child_ids,i)
    fi = getindex!(fi_cache,fields,id)
    mi = getindex!(mi_cache,cmaps,id)
    zi = Fields.evaluate!(zi_cache,mi,xi)
    y_cache.array[i] = Fields.evaluate!(yi_cache,fi,zi)
  end
  return y_cache.array
end

"""
function Geometry.return_cache(a::FineToCoarseField,x::Table{<:Point})
  @check length(x) == num_subcells(a.rrule)
  fields, rr = a.fine_fields, a.rrule
  cmaps = get_inverse_cell_map(a.rrule)

  xi_cache = array_cache(x.data)
  fi_cache = array_cache(fields)
  mi_cache = array_cache(cmaps)

  xi = getindex!(xi_cache,x.data,1)
  child_id = x_to_cell(rr,xi)
  mi = getindex!(mi_cache,cmaps,child_id)
  fi = getindex!(fi_cache,fields,child_id)

  zi_cache = Fields.return_cache(mi,xi)
  zi = evaluate!(zi_cache,mi,xi)

  yi_type  = Fields.return_type(fi,zi)
  yi_cache = Fields.return_cache(fi,zi)
  y_cache  = Arrays.CachedArray(yi_type,1)

  return fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache
end

function Geometry.evaluate!(cache,a::FineToCoarseField,x::Table{<:Point})
  @check length(x) == num_subcells(a.rrule)
  fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache = cache
  cmaps  = get_inverse_cell_map(a.rrule)
  fields = a.fine_fields

  Arrays.setsize!(y_cache, size(x.data))
  for i in eachindex(x)
    for j in x.ptrs[i]:x.ptrs[i+1]
      xi = getindex!(xi_cache,x.data,j)
      fi = getindex!(fi_cache,fields,i)
      mi = getindex!(mi_cache,cmaps,i)
      zi = Fields.evaluate!(zi_cache,mi,xi)
      y_cache.array[j] = Fields.evaluate!(yi_cache,fi,zi)
    end
  end
  return y_cache.array
end
"""