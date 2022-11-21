
struct RefinementRule{T,A,B,C,D}
  reffe::T
  ref_grid::A
  f2c_cell_map::B
  c2f_cell_map::C
  measures::D
end

function RefinementRule(reffe::LagrangianRefFE{D},nrefs::Integer; kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(reffe,partition;kwargs...)
end

function RefinementRule(reffe::LagrangianRefFE{D},partition::NTuple{D,Integer}; kwargs...) where D
  ref_grid = UnstructuredGrid(compute_reference_grid(reffe,partition))
  f2c_cell_map = Geometry.get_cell_map(ref_grid)
  c2f_cell_map = lazy_map(Fields.inverse_map,f2c_cell_map)

  ref_trian = Triangulation(UnstructuredDiscreteModel(ref_grid))
  measures  = get_cell_measure(ref_trian)
  M = sum(measures)
  measures /= M

  return RefinementRule(reffe,ref_grid,f2c_cell_map,c2f_cell_map,measures)
end

ReferenceFEs.get_polytope(rr::RefinementRule) = ReferenceFEs.get_polytope(rr.reffe)
get_ref_grid(rr::RefinementRule) = rr.ref_grid
get_measures(rr::RefinementRule) = rr.measures

function Geometry.get_cell_map(rr::RefinementRule)
  return rr.f2c_cell_map
end

function get_inverse_cell_map(rr::RefinementRule)
  return rr.c2f_cell_map
end
