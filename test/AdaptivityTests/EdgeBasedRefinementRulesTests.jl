module EdgeBasedRefinementRulesTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays

D = 2
sol(x) = x[1] + x[2]

rr = Gridap.Adaptivity.GreenRefinementRule(QUAD,3)
rr_model = Gridap.Adaptivity.get_ref_grid(rr)

# Checking that inv_cmaps ∘ cmaps = identity
cmaps = get_cell_map(rr_model)
inv_cmaps = lazy_map(Gridap.Fields.inverse_map,cmaps)

pts = map(x -> VectorValue(rand(D)),1:10)
for p in pts
  ichild = Gridap.Adaptivity.x_to_cell(rr,p)
  m = cmaps[ichild]
  m_inv = inv_cmaps[ichild]

  y = evaluate(m,p)
  z = evaluate(m_inv,y)
  println(ichild, " :: ", p," -> ",y, " -> ", z, " - ", p ≈ z)
end

# This is the problematic cellmap
cm1 = cmaps[1]

# Model topology: The coordinates for each cell seem to be fine....
cell_node_ids = rr_model.grid.cell_node_ids
node_coords   = rr_model.grid.node_coordinates
polytopes     = get_polytopes(rr_model)
cell_type     = get_cell_type(rr_model)
orientation   = OrientationStyle(rr_model)

cell_coords   = map(ids->node_coords[ids],cell_node_ids)

# Even if we re-build the UnstructuredDiscreteModel from scratch, the problem persists....
_topo   = UnstructuredGridTopology(node_coords,cell_node_ids,cell_type,polytopes,orientation)
_reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(_topo))
_grid   = UnstructuredGrid(get_vertex_coordinates(_topo),get_faces(_topo,D,0),_reffes,get_cell_type(_topo),OrientationStyle(_topo))
_labels = FaceLabeling(_topo)
_model  = UnstructuredDiscreteModel(_grid,_topo,_labels)

_cmaps = get_cell_map(_model)
_inv_cmaps = lazy_map(Gridap.Fields.inverse_map,cmaps)

_trian = Triangulation(_model)
_p2c_cache = CellData._point_to_cell_cache(CellData.KDTreeSearch(),_trian)
x_to_cell(x) = CellData._point_to_cell!(_p2c_cache,x)

pts = map(x -> VectorValue(rand(D)),1:10)
for p in pts
  ichild = x_to_cell(p)
  m = cmaps[ichild]
  m_inv = inv_cmaps[ichild]

  y = evaluate(m,p)
  z = evaluate(m_inv,y)
  println(ichild, " :: ", p," -> ",y, " -> ", z, " - ", p ≈ z)
end