module TestVisualizationData

using Test
using Gridap.Geometry
using Gridap.CellData
using Gridap.TensorValues
using Gridap.Algebra
using Gridap.Visualization
using Gridap.ReferenceFEs
using Gridap.Geometry: get_reffes

domain = (0,2pi, 0, 4pi)
partition = (2, 3)
model = simplexify(CartesianDiscreteModel(domain,partition))
f(pt) = sin(pt[1])

trian = Triangulation(model)
u = CellField(f,trian)

visdata = visualization_data(trian, cellfields=Dict("u" =>u,"f"=>f))

ncells = num_cells(visdata.grid)
nodes_per_cell = num_nodes(first(get_reffes(visdata.grid)))
@test size(visdata.nodaldata["u"]) == (ncells * nodes_per_cell,)
@test size(visdata.celldata["cell"]) == (ncells,)
@test visdata.grid isa Visualization.VisualizationGrid
@test visdata isa Visualization.VisualizationData

end#module
