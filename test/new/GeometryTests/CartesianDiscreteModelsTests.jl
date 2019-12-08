module CartesianDiscreteModelsTests

using Gridap
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: CartesianMap
using Test

domain = (0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)
test_discrete_model(model)

map = get_cell_map(model)
@test isa(map,CartesianMap)

@test get_polytope(model) == HEX
@test get_polytope(Polytope{1},model) == SEGMENT

function polar(q)
  r, t, z = q
  x = r*cos(t)
  y = r*sin(t)
  Point(x,y,z)
end

domain = (1,2,0,pi,0,0.5)
partition = (10,30,4)
model = CartesianDiscreteModel(domain,partition,polar)
test_discrete_model(model)

#using Gridap.Visualization
#writevtk(model,"model")

end # module
