module InverseCellMapsTests

using Gridap
using Gridap.Arrays
using Gridap.CellData
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs

# Build model
node_coords = VectorValue{2, Float64}[
  VectorValue{2, Float64}(0.0, 0.0),
  VectorValue{2, Float64}(1.0, 0.0),
  VectorValue{2, Float64}(0.0, 1.0),
  VectorValue{2, Float64}(1.0, 1.0),
  VectorValue{2, Float64}(0.0, 0.5)
]

data = Int32[2, 4, 5,
             1, 2, 5,
             3, 4, 5]
ptrs = Int32[1,4,7,10]
cell_node_ids = Table(data,ptrs)

polytopes = [TRI]
cell_type = Int8[1,1,1]
orientation = Oriented()

topo   = UnstructuredGridTopology(node_coords,cell_node_ids,cell_type,polytopes,orientation)
reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
grid   = UnstructuredGrid(get_vertex_coordinates(topo),get_faces(topo,2,0),reffes,get_cell_type(topo),OrientationStyle(topo))
labels = FaceLabeling(topo)
model  = UnstructuredDiscreteModel(grid,topo,labels)

cmaps = get_cell_map(model)
inv_cmaps = lazy_map(Gridap.Fields.inverse_map,cmaps)

trian = Triangulation(model)
p2c_cache = CellData._point_to_cell_cache(CellData.KDTreeSearch(),trian)
x_to_cell(x) = CellData._point_to_cell!(p2c_cache,x)

"""
 As can be seen in the loop, the inverse map fails for the first cell of the model. 
"""
pts = map(x -> VectorValue(rand(2)),1:10)
for p in pts
  ichild = x_to_cell(p)
  m = cmaps[ichild]
  m_inv = inv_cmaps[ichild]

  y = evaluate(m_inv,p) # QUAD -> Subelement
  z = evaluate(m,y)     # Subelement -> QUAD
  println(ichild, " :: ", p," -> ",y, " -> ", z, " - ", p ≈ z)
end

end