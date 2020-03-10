module QPointCellFields

using Gridap.Geometry
using Gridap.Geometry: ArrayOfEvaluatedFields
using Gridap.Fields
using Gridap.Integration

domain = (0,1,0,1)
n = 3
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

trian = Triangulation(model)
degree = 4
quad = CellQuadrature(trian,degree)

q = get_coordinates(quad)
s_q = [i*ones(size(qi)) for (i,qi) in enumerate(q)]
a = ArrayOfEvaluatedFields(s_q,q)
test_array_of_fields(a,q,s_q)

r = QPointCellField(0.0,trian,quad)
r_q = evaluate(r,q)
test_cell_field(r,q,r_q)

end # module
