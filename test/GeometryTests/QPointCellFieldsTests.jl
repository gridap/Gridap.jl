module QPointCellFields

using Gridap.Geometry
using Gridap.Fields
using Gridap.Integration
using Gridap.CellData

domain = (0,1,0,1)
n = 3
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

trian = Triangulation(model)
degree = 4
quad = CellQuadrature(trian,degree)
x = get_coordinates(quad)

r = CellField(0.0,quad)
r_x = evaluate(r,x)
test_cell_field(r,x,r_x)

end # module
