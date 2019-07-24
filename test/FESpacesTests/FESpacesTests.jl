module FESpacesTests

using Test
using Gridap

model = CartesianDiscreteModel(partition=(2,2))

order = 1
tags = [1,2,3,4]
fespace = ConformingFESpace(Float64,model,order,tags)

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

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

end # module

