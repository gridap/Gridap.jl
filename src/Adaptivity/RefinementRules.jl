
struct RefinementRule{T,A,B}
  reffe::T
  f2c_cell_map::A
  c2f_cell_map::B
end

function RefinementRule(reffe::LagrangianRefFE{D},nrefs::Integer; kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(reffe,partition;kwargs...)
end

function RefinementRule(reffe::LagrangianRefFE{D},partition::NTuple{D,Integer}; kwargs...) where D
  ref_grid = UnstructuredGrid(compute_reference_grid(reffe,partition))
  f2c_cell_map = Geometry.get_cell_map(ref_grid)
  c2f_cell_map = lazy_map(Fields.inverse_map,f2c_cell_map)
  return RefinementRule(reffe,f2c_cell_map,c2f_cell_map)
end

function Geometry.get_cell_map(rr::RefinementRule)
  return rr.f2c_cell_map
end

function get_inverse_cell_map(rr::RefinementRule)
  return rr.c2f_cell_map
end
