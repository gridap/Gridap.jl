
abstract type RefinementRuleType end
struct GenericRefinement <: RefinementRuleType end

struct RefinementRule{P,T,A,B,C}
  poly         :: P
  ref_grid     :: A
  f2c_cell_map :: B
  x_to_cell    :: C

  function RefinementRule{T}(poly         ::Polytope,
                             ref_grid     ::Grid,
                             f2c_cell_map ::AbstractVector,
                             x_to_cell    ::Function) where T <: RefinementRuleType
    P = typeof(poly)
    A = typeof(ref_grid)
    B = typeof(f2c_cell_map)
    C = typeof(x_to_cell)
    new{P,T,A,B,C}(poly,ref_grid,f2c_cell_map,x_to_cell)
  end
end

ReferenceFEs.get_polytope(rr::RefinementRule) = rr.poly
get_ref_grid(rr::RefinementRule) = rr.ref_grid
num_subcells(rr::RefinementRule) = num_cells(rr.ref_grid)

function Geometry.get_cell_map(rr::RefinementRule)
  return rr.f2c_cell_map
end

function get_inverse_cell_map(rr::RefinementRule)
  return lazy_map(Fields.inverse_map,rr.f2c_cell_map)
end

function get_cell_measures(rr::RefinementRule)
  ref_grid  = get_ref_grid(rr)
  ref_trian = Triangulation(UnstructuredDiscreteModel(ref_grid))
  measures  = get_cell_measure(ref_trian)
  M = sum(measures)
  measures /= M
  return measures
end

function RefinementRule(::T,poly::Polytope, ref_grid::Grid; cell_maps=nothing, x_to_cell_id=nothing) where T <: RefinementRuleType
  if (cell_maps === nothing)
    f2c_cell_map = Geometry.get_cell_map(ref_grid)
  else
    @check length(cell_maps) == num_cells(ref_grid)
    f2c_cell_map = cell_maps
  end

  if x_to_cell_id === nothing
    ref_trian    = Triangulation(UnstructuredDiscreteModel(ref_grid))
    p2c_cache    = CellData._point_to_cell_cache(CellData.KDTreeSearch(),ref_trian)
    x_to_cell(x) = CellData._point_to_cell!(p2c_cache, x)
  else
    x_to_cell = x_to_cell_id
  end

  return RefinementRule{T}(poly,ref_grid,f2c_cell_map,x_to_cell)
end

# GenericRefinement Rule

function RefinementRule(reffe::LagrangianRefFE{D},nrefs::Integer;kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(get_polytope(reffe),partition;kwargs...)
end

function RefinementRule(reffe::LagrangianRefFE{D},partition::NTuple{D,Integer};kwargs...) where D
  return RefinementRule(get_polytope(reffe),partition;cell_maps=cell_maps)
end

function RefinementRule(poly::Polytope{D},nrefs::Integer;kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(poly,partition;kwargs...)
end

function RefinementRule(poly::Polytope{D},partition::NTuple{D,Integer};kwargs...) where D
  ref_grid = UnstructuredGrid(compute_reference_grid(poly,partition))
  return RefinementRule(GenericRefinement(),poly,ref_grid;kwargs...)
end

# FineToCoarseFields

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
  fields, cmaps = a.fine_fields, a.rrule.c2f_cell_map

  fi_cache = array_cache(fields)
  cm_cache = array_cache(cmaps)
  xi_cache = Fields.return_cache(first(cmaps),x)
  yi_cache = Fields.return_cache(first(fields),x)
  return fi_cache, cm_cache, xi_cache, yi_cache
end

function Geometry.evaluate!(cache,a::FineToCoarseField,x::Point)
  fi_cache, cm_cache, xi_cache, yi_cache = cache
  fields, x_to_cell, cmaps = a.fine_fields, a.rrule.x_to_cell, a.rrule.c2f_cell_map

  child_id = x_to_cell(x) # Find correct subcell
  fi = getindex!(fi_cache,fields,child_id)
  mi = getindex!(cm_cache,cmaps,child_id)
  xi = Fields.evaluate!(xi_cache,mi,x)  # xc -> xf
  yi = Fields.evaluate!(yi_cache,fi,xi) # xf -> yf
  return yi
end

function Geometry.return_cache(a::FineToCoarseField,x::AbstractArray{<:Point})
  fields, x_to_cell, cmaps = a.fine_fields, a.rrule.x_to_cell, a.rrule.c2f_cell_map

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
  fields, x_to_cell, cmaps = a.fine_fields, a.rrule.x_to_cell, a.rrule.c2f_cell_map

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
