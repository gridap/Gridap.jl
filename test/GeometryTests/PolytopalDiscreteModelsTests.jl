module PolytopalDiscreteModelsTests
using Test
using Gridap
using Gridap.Geometry, Gridap.ReferenceFEs, Gridap.TensorValues
using Gridap.Fields, Gridap.FESpaces, Gridap.CellData
using Gridap.Io

using FillArrays

function test_model(model)
  D = num_cell_dims(model)

  # Model API
  test_discrete_model(model)
  Geometry.restrict(model,[1,2,3,4])

  # Triangulations

  for d in 0:D
    trian = Triangulation(ReferenceFE{d},model)
    test_triangulation(trian)
  end

  Ω = Triangulation(ReferenceFE{D},model,[1,2,3,4])
  Γ = Boundary(model)
  Λ = Skeleton(model)

  d = tempdir()
  writevtk(model,d*"/polytopal_model")
  writevtk(Ω,d*"/polytopal_trian")
  writevtk(Γ,d*"/polytopal_boundary")
  writevtk(Λ,d*"/polytopal_skeleton")
end

model = CartesianDiscreteModel((0,1,0,1),(2,2))

pmodel = Geometry.PolytopalDiscreteModel(model)
pgrid = Geometry.PolytopalGrid(get_polytopes(pmodel))
test_model(pmodel)

vmodel = Geometry.voronoi(Geometry.simplexify(model))
test_model(vmodel)

model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
pmodel = Geometry.PolytopalDiscreteModel(model)
test_model(pmodel)

# IO

model = CartesianDiscreteModel((0,1,0,1),(2,2))
topo = Geometry.PolytopalGridTopology(get_grid_topology(model))

dict = to_dict(topo)
@test dict[:Dc] == 2
@test dict[:Dp] == 2
@test haskey(dict,:cell_vertices)
@test haskey(dict,:polytopes)
topo_r = from_dict(Geometry.PolytopalGridTopology,dict)
@test num_faces(topo_r,0) == num_faces(topo,0)
@test num_faces(topo_r,2) == num_faces(topo,2)
@test get_vertex_coordinates(topo_r) ≈ get_vertex_coordinates(topo)
cell_verts_orig = get_faces(topo,2,0)
cell_verts_rest = get_faces(topo_r,2,0)
@test cell_verts_orig.data == cell_verts_rest.data
@test cell_verts_orig.ptrs == cell_verts_rest.ptrs

grid = Geometry.PolytopalGrid(get_grid(model))
dict = to_dict(grid)
@test dict[:Dc] == 2
@test dict[:Dp] == 2
grid_r = from_dict(Geometry.PolytopalGrid,dict)
@test num_cells(grid_r) == num_cells(grid)
@test num_nodes(grid_r) == num_nodes(grid)
@test get_node_coordinates(grid_r) ≈ get_node_coordinates(grid)
ids_orig = get_cell_node_ids(grid)
ids_rest = get_cell_node_ids(grid_r)
@test ids_orig.data == ids_rest.data
@test ids_orig.ptrs == ids_rest.ptrs

model = CartesianDiscreteModel((0,1,0,1),(3,3))
pmodel = Geometry.PolytopalDiscreteModel(model)
dict = to_dict(pmodel)
@test haskey(dict,:grid)
@test haskey(dict,:topology)
@test haskey(dict,:labeling)
pmodel_r = from_dict(Geometry.PolytopalDiscreteModel,dict)
@test num_cells(pmodel_r) == num_cells(pmodel)
@test num_faces(get_grid_topology(pmodel_r),0) == num_faces(get_grid_topology(pmodel),0)
@test num_faces(get_grid_topology(pmodel_r),1) == num_faces(get_grid_topology(pmodel),1)
@test num_faces(get_grid_topology(pmodel_r),2) == num_faces(get_grid_topology(pmodel),2)
@test get_node_coordinates(get_grid(pmodel_r)) ≈ get_node_coordinates(get_grid(pmodel))

vmodel = Geometry.voronoi(model)
dict_v = to_dict(vmodel)
vmodel_r = from_dict(Geometry.PolytopalDiscreteModel,dict_v)
@test num_cells(vmodel_r) == num_cells(vmodel)
@test num_faces(get_grid_topology(vmodel_r),0) == num_faces(get_grid_topology(vmodel),0)

model = CartesianDiscreteModel((0,1,0,1),(2,2))
pmodel = Geometry.PolytopalDiscreteModel(model)
json_str = to_json(pmodel)
pmodel_rj = from_json(Geometry.PolytopalDiscreteModel,json_str)
@test num_cells(pmodel_rj) == num_cells(pmodel)

end