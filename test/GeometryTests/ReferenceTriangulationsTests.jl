module RestrictedTriangulationsTests

using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Test

domain = (-1,1,-1,1)
partition = (10,10)
model = CartesianDiscreteModel(domain,partition)
trian = get_triangulation(model)

test_triangulation(ReferenceTriangulation(trian))


end # module
