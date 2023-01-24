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

cell_node_ids = rr_model.grid.cell_node_ids
node_coords = rr_model.grid.node_coordinates
cell_coords = map(ids->node_coords[ids],cell_node_ids)

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

cm1 = cmaps[1]