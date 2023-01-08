
abstract type RefinementRuleType end
struct GenericRefinement <: RefinementRuleType end
struct WithoutRefinement <: RefinementRuleType end

struct RefinementRule{T,P,A,B,C}
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
    new{T,P,A,B,C}(poly,ref_grid,f2c_cell_map,x_to_cell)
  end
end

ReferenceFEs.get_polytope(rr::RefinementRule) = rr.poly
get_ref_grid(rr::RefinementRule) = rr.ref_grid
num_subcells(rr::RefinementRule) = num_cells(rr.ref_grid)
num_ref_faces(rr::RefinementRule,d::Int) = num_faces(rr.ref_grid,d)
RefinementRuleType(::RefinementRule{T}) where T = T

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

function bundle_points_by_subcell(rr::RefinementRule,x::AbstractArray{<:Point})
  x_to_cell = rr.x_to_cell
  npts      = length(x)
  nchildren = num_subcells(rr)

  child_ids = map(x_to_cell,x)
  ptrs      = fill(0,nchildren+1)
  for i in 1:npts
    ptrs[child_ids[i]+1] += 1
  end
  ptrs[1] = 1
  
  data = lazy_map(Reindex(x),sortperm(child_ids))
  return Table(data,ptrs)
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

