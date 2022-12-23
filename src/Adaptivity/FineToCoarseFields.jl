
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
  fields, x_to_cell = a.fine_fields, a.rrule.x_to_cell
  cmaps = get_inverse_cell_map(a.rrule)

  child_id = x_to_cell(x) # Find correct subcell
  fi = getindex!(fi_cache,fields,child_id)
  mi = getindex!(cm_cache,cmaps,child_id)
  xi = Fields.evaluate!(xi_cache,mi,x)  # xc -> xf
  yi = Fields.evaluate!(yi_cache,fi,xi) # xf -> yf
  return yi
end

function Geometry.return_cache(a::FineToCoarseField,x::AbstractArray{<:Point})
  fields, x_to_cell = a.fine_fields, a.rrule.x_to_cell
  cmaps = get_inverse_cell_map(a.rrule)

  xi_cache = array_cache(x)
  fi_cache = array_cache(fields)
  mi_cache = array_cache(cmaps)

  xi = getindex!(xi_cache,x,1)
  child_id = x_to_cell(xi)
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
  fields, x_to_cell = a.fine_fields, a.rrule.x_to_cell
  cmaps = get_inverse_cell_map(a.rrule)

  Arrays.setsize!(y_cache, size(x))

  for i in eachindex(x)
    xi = getindex!(xi_cache,x,i)
    child_id = x_to_cell(xi)
    fi = getindex!(fi_cache,fields,child_id)
    mi = getindex!(mi_cache,cmaps,child_id)
    zi = Fields.evaluate!(zi_cache,mi,xi)
    y_cache.array[i] = Fields.evaluate!(yi_cache,fi,zi)
  end
  return y_cache.array
end
