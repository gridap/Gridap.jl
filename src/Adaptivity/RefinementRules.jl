
abstract type RefinementRuleType end
struct GenericRefinement <: RefinementRuleType end
struct WithoutRefinement <: RefinementRuleType end

struct RefinementRule{P,A}
  T            :: RefinementRuleType
  poly         :: P
  ref_grid     :: A
  p2c_cache    :: Tuple
end

function RefinementRule(T::RefinementRuleType,poly::Polytope,ref_grid::Grid)
  ref_trian = Triangulation(UnstructuredDiscreteModel(ref_grid))
  p2c_cache = CellData._point_to_cell_cache(CellData.KDTreeSearch(),ref_trian)
  return RefinementRule(T,poly,ref_grid,p2c_cache)
end

function Base.show(io::IO,rr::RefinementRule{P}) where P
  T = RefinementRuleType(rr)
  print(io,"RefinementRule{$T,$P}")
end

ReferenceFEs.get_polytope(rr::RefinementRule) = rr.poly
get_ref_grid(rr::RefinementRule) = rr.ref_grid
num_subcells(rr::RefinementRule) = num_cells(rr.ref_grid)
num_ref_faces(rr::RefinementRule,d::Int) = num_faces(rr.ref_grid,d)
RefinementRuleType(rr::RefinementRule) :: RefinementRuleType = rr.T

function Geometry.get_cell_map(rr::RefinementRule)
  ref_grid = get_ref_grid(rr)
  return Geometry.get_cell_map(ref_grid)
end

function get_inverse_cell_map(rr::RefinementRule)
  f2c_cell_map = get_cell_map(rr)
  return lazy_map(Fields.inverse_map,f2c_cell_map)
end

function get_cell_measures(rr::RefinementRule)
  ref_grid  = get_ref_grid(rr)
  ref_trian = Triangulation(UnstructuredDiscreteModel(ref_grid))
  measures  = get_cell_measure(ref_trian)
  M = sum(measures)
  measures /= M
  return measures
end

function get_cell_polytopes(rr::RefinementRule)
  ref_grid   = get_ref_grid(rr)
  polys      = get_polytopes(ref_grid)
  cell_types = get_cell_type(ref_grid)
  return CompressedArray(polys,cell_types)
end

function get_cell_polytopes(rrules::AbstractArray{<:RefinementRule})
  rr_polys = lazy_map(rr->get_polytopes(rr.ref_grid),rrules)

  # NOTE: The innermost `unique` is to optimize for CompressedArrays
  polys_new = unique(reduce(vcat,unique(rr_polys)))

  rr_cell_type = lazy_map(rr->get_cell_type(rr.ref_grid),rrules)
  rr2new_cell_type  = lazy_map(vp->map(p->findfirst(x->x==p,polys_new),vp),rr_polys)
  cell_type_new = reduce(vcat,lazy_map((gids,lids)->lazy_map(Reindex(gids),lids),rr2new_cell_type,rr_cell_type))

  return polys_new, cell_type_new
end

x_to_cell(rr::RefinementRule,x::Point) = CellData._point_to_cell!(rr.p2c_cache,x)

function bundle_points_by_subcell(rr::RefinementRule,x::AbstractArray{<:Point})
  npts      = length(x)
  nchildren = num_subcells(rr)

  child_ids = map(xi -> x_to_cell(rr,xi),x)
  ptrs      = fill(0,nchildren+1)
  for i in 1:npts
    ptrs[child_ids[i]+1] += 1
  end
  ptrs[1] = 1
  
  data = lazy_map(Reindex(x),sortperm(child_ids))
  return Table(data,ptrs)
end

# GenericRefinement Rule

function RefinementRule(reffe::LagrangianRefFE{D},nrefs::Integer;kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(get_polytope(reffe),partition;kwargs...)
end

function RefinementRule(reffe::LagrangianRefFE{D},partition::NTuple{D,Integer};kwargs...) where D
  return RefinementRule(get_polytope(reffe),partition;kwargs...)
end

function RefinementRule(poly::Polytope{D},nrefs::Integer;kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(poly,partition;kwargs...)
end

function RefinementRule(poly::Polytope{D},partition::NTuple{D,Integer};kwargs...) where D
  ref_grid = UnstructuredGrid(compute_reference_grid(poly,partition))
  return RefinementRule(GenericRefinement(),poly,ref_grid;kwargs...)
end

