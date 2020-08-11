module TestVisualizationData

using Test
using Gridap
using Gridap.Visualization
using Gridap.ReferenceFEs
using Gridap.Geometry: get_reffes

domain = (0,2pi, 0, 4pi)
partition = (2, 3)
model = simplexify(CartesianDiscreteModel(domain,partition))
f(pt) = sin(pt[1])
V0 = TestFESpace(
  reffe=:Lagrangian, order=1, valuetype=Float64,
  conformity=:H1, model=model)
U = TrialFESpace(V0)
trian = Triangulation(model)
degree = 2
quad = CellQuadrature(trian,degree)
A(u, v) = u ⊙ v
b(v) = v⊙f
t_Ω = AffineFETerm(A,b,trian,quad)
op = AffineFEOperator(U,V0,t_Ω)
u = solve(op)

visdata = visualization_data(trian, cellfields=Dict("u" =>u))

ncells = num_cells(visdata.grid)
nodes_per_cell = num_nodes(first(get_reffes(visdata.grid)))
@test size(visdata.nodaldata["u"]) == (ncells * nodes_per_cell,)
@test size(visdata.celldata["cell"]) == (ncells,)
@test visdata.grid isa Gridap.Visualization.VisualizationGrid
@test visdata isa Visualization.VisualizationData

end#module
