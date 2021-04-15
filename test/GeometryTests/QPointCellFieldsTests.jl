module QPointCellFields

using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData

domain = (0,1,0,1)
n = 3
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

trian = Triangulation(model)
degree = 4
quad = CellQuadrature(trian,degree)
q = get_coordinates(quad)

r = CellField(0.0,trian,quad)
r_q = evaluate(r,q)
test_cell_field(r,q,r_q)

end # module
