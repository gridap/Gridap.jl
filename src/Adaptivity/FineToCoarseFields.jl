
"""
    struct FineToCoarseField <: Field
      fine_fields :: AbstractVector{<:Field}
      rrule       :: RefinementRule
      id_map      :: AbstractVector{<:Integer}
    end

Given a domain and a non-overlapping refined cover, a `FineToCoarseField`
is a `Field` defined in the domain and constructed by a set of fields defined on 
the subparts of the covering partition. The refined cover is represented by a `RefinementRule`. 

# Parameters: 

- `rrule`: Refinement rule representing the covering partition.
- `fine_fields`: Fields defined on the subcells of the covering partition. To accomodate the
                 case where not all subcells are present (e.g in distributed), we allow for 
                 `length(fine_fields) != num_subcells(rrule)`.
- `id_map`: Mapping from the subcell ids to the indices of the `fine_fields` array. 
            If `length(fine_fields) == num_subcells(rrule)`, this is the identity map. 
            Otherwise, the `id_map` is used to map the subcell ids to the indices of the 
            `fine_fields` array.
"""
struct FineToCoarseField{A,B,C} <: Field
  fine_fields :: A
  rrule       :: B
  id_map      :: C
  function FineToCoarseField(
    fine_fields::AbstractArray{<:Field},
    rrule::RefinementRule,
    child_ids::AbstractArray{<:Integer}
  )
    n_children = num_subcells(rrule)
    id_map = Arrays.find_inverse_index_map(child_ids,n_children)
    A = typeof(fine_fields)
    B = typeof(rrule)
    C = typeof(id_map)
    new{A,B,C}(fine_fields,rrule,id_map)
  end
  function FineToCoarseField(
    fine_fields::AbstractArray{<:Field},
    rrule::RefinementRule,
  )
    n_children = num_subcells(rrule)
    id_map = Base.OneTo(n_children)
    @check length(fine_fields) == n_children
    A = typeof(fine_fields)
    B = typeof(rrule)
    C = typeof(id_map)
    new{A,B,C}(fine_fields,rrule,id_map)
  end
end

Arrays.return_type(a::FineToCoarseField,x::Point) = return_type(first(a.fine_fields),zero(x))
Arrays.return_value(a::FineToCoarseField,x::Point) = zero(return_type(a,x))
Arrays.return_cache(a::FineToCoarseField,x::Point) = return_cache(a,[x])
Arrays.evaluate!(cache,a::FineToCoarseField,x::Point) = first(evaluate!(cache,a,[x]))

function Arrays.return_type(a::FineToCoarseField,x::AbstractArray{<:Point})
  T = return_type(first(a.fine_fields),zero(eltype(x)))
  return Vector{T}
end

function Arrays.return_value(a::FineToCoarseField,x::AbstractArray{<:Point})
  T = return_type(a,x)
  return similar(T,size(x))
end

function Arrays.return_cache(a::FineToCoarseField,x::AbstractArray{<:Point})
  fields, rr, id_map = a.fine_fields, a.rrule, a.id_map
  cmaps = get_inverse_cell_map(rr)

  xi_cache = array_cache(x)
  fi_cache = array_cache(fields)
  mi_cache = array_cache(cmaps)

  # Geometric caches
  xi = first(x)
  mi = getindex!(mi_cache,cmaps,1)
  zi_cache = Fields.return_cache(mi,xi)
  zi = zero(Fields.return_type(mi,xi))

  # Evaluation caches
  fi = getindex!(fi_cache,fields,1)
  yi_cache = Fields.return_cache(fi,zi)
  
  # Output cache
  yi_types = map(fii -> Fields.return_type(fii,zi), fields)
  yi_type  = first(yi_types)
  @assert all(t -> yi_type == t, yi_types)
  y_cache  = Arrays.CachedArray(zeros(yi_type,size(x)))

  return fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache
end

function Arrays.evaluate!(cache,a::FineToCoarseField,x::AbstractArray{<:Point})
  fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache = cache
  fields, rr, id_map = a.fine_fields, a.rrule, a.id_map
  cmaps = get_inverse_cell_map(rr)

  Arrays.setsize!(y_cache, size(x))
  y = y_cache.array
  fill!(y,zero(eltype(y)))

  for i in eachindex(x)
    xi = getindex!(xi_cache,x,i)
    child_id = x_to_cell(rr,xi)
    field_id = id_map[child_id]
    if field_id != 0
      fi = getindex!(fi_cache,fields,field_id)
      mi = getindex!(mi_cache,cmaps,child_id)
      zi = Fields.evaluate!(zi_cache,mi,xi)
      y[i] = Fields.evaluate!(yi_cache,fi,zi)
    end
  end
  return y
end

# Fast evaluation of FineToCoarseFields: 
# Points are pre-classified into the children cells, which allows for the search to be 
# skipped entirely.

function Arrays.return_type(a::FineToCoarseField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  return_type(a,x)
end

function Arrays.return_value(a::FineToCoarseField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  return_value(a,x)
end

function Arrays.return_cache(a::FineToCoarseField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  return_cache(a,x)
end

function Arrays.evaluate!(cache,a::FineToCoarseField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache = cache
  fields, rr, id_map = a.fine_fields, a.rrule, a.id_map
  cmaps = get_inverse_cell_map(rr)

  Arrays.setsize!(y_cache, size(x))
  y = y_cache.array
  fill!(y,zero(eltype(y)))

  for i in eachindex(x)
    xi = getindex!(xi_cache,x,i)
    child_id = child_ids[i]
    field_id = id_map[child_id]
    if field_id != 0
      fi = getindex!(fi_cache,fields,field_id)
      mi = getindex!(mi_cache,cmaps,field_id)
      zi = Fields.evaluate!(zi_cache,mi,xi)
      y[i] = Fields.evaluate!(yi_cache,fi,zi)
    end
  end
  return y
end
