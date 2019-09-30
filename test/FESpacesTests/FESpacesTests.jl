module FESpacesTests

using Test
using Gridap

model = CartesianDiscreteModel(partition=(2,2))

order = 1
tags = [1,2,3,4]
fespace = H1ConformingFESpace(Float64,model,order,tags)

@test string(fespace) == "ConformingFESpace object"
s = "ConformingFESpace object:\n physdim: 2\n refdim: 2\n valuetype: Float64"
@test sprint(show,"text/plain",fespace) == s

trian = Triangulation(model)
quad = CellQuadrature(trian,degree=2)

cellbasis = CellBasis(fespace)
a(v,u) = varinner(v,u)
bfun(x) = x[2]
b(v) = varinner(v,CellField(trian,bfun))
mmat = integrate(a(cellbasis,cellbasis),trian,quad)
bvec = integrate(b(cellbasis),trian,quad)

fun(x) = x[1] + x[2]

test_fe_space(fespace,5,4,mmat,bvec,fun)

V = TestFESpace(fespace)
test_fe_space(V,5,4,mmat,bvec,fun)

U = TrialFESpace(fespace,fun)
test_fe_space(U,5,4,mmat,bvec,fun)

fun(x) = x[1] + x[2]
z = 2.0
U = TrialFESpace(fespace,z)
test_fe_space(U,5,4,mmat,bvec,fun)

U = TrialFESpace(fespace,[z,z,z,z])
test_fe_space(U,5,4,mmat,bvec,fun)

end # module
