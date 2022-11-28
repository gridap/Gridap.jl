
struct RefinementRule{T,A,B,C,D,E}
  poly         :: T
  ref_grid     :: A
  f2c_cell_map :: B
  c2f_cell_map :: C
  measures     :: D
  x_to_cell    :: E
end

function RefinementRule(reffe::LagrangianRefFE{D},nrefs::Integer;kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(reffe,partition;kwargs...)
end

function RefinementRule(reffe::LagrangianRefFE{D},partition::NTuple{D,Integer}) where D
  ref_grid  = UnstructuredGrid(compute_reference_grid(reffe,partition))
  cell_maps = Geometry.get_cell_map(ref_grid)
  return RefinementRule(get_polytope(reffe),partition;cell_maps=cell_maps)
end

function RefinementRule(poly::Polytope{D},nrefs::Integer;kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(poly,partition;kwargs...)
end

function RefinementRule(poly::Polytope{D},partition::NTuple{D,Integer};cell_maps=nothing) where D
  ref_grid = UnstructuredGrid(compute_reference_grid(poly,partition))
  if (cell_maps === nothing)
    f2c_cell_map = Geometry.get_cell_map(ref_grid)
  else
    @check length(cell_maps) == num_cells(ref_grid)
    f2c_cell_map = cell_maps
  end
  c2f_cell_map = lazy_map(Fields.inverse_map,f2c_cell_map)

  ref_trian = Triangulation(UnstructuredDiscreteModel(ref_grid))
  measures  = get_cell_measure(ref_trian)
  M = sum(measures)
  measures /= M

  p2c_cache    = CellData._point_to_cell_cache(CellData.KDTreeSearch(),ref_trian)
  x_to_cell(x) = CellData._point_to_cell!(p2c_cache, x)

  return RefinementRule(poly,ref_grid,f2c_cell_map,c2f_cell_map,measures,x_to_cell)
end

ReferenceFEs.get_polytope(rr::RefinementRule) = rr.poly
get_ref_grid(rr::RefinementRule) = rr.ref_grid
get_measures(rr::RefinementRule) = rr.measures
num_subcells(rr::RefinementRule) = length(rr.measures)

function Geometry.get_cell_map(rr::RefinementRule)
  return rr.f2c_cell_map
end

function get_inverse_cell_map(rr::RefinementRule)
  return rr.c2f_cell_map
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
