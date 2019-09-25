module ZeroMeanFESpacesTests

using Test
using Gridap

model = CartesianDiscreteModel(partition=(2,2))

order = 2

fespace = FESpace(
  reffe = :PLagrangian,
  conformity = :L2,
  valuetype = Float64,
  order = order,
  model = model)

quadorder = order
fespace = ZeroMeanFESpace(fespace,quadorder)

@test fespace.vol ≈ 4.0

ufun(x) = x[1] + x[2] + 2.0

trian = Triangulation(model)
quad = CellQuadrature(trian,degree=order)

bh = FEBasis(fespace)
uh = zero(fespace)

cellmat = integrate(inner(bh,bh),trian,quad)
cellvec = integrate(inner(bh,uh),trian,quad)

nfree = 23

ndiri = 0

test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

V = TrialFESpace(fespace)

x = rand(Float64,num_free_dofs(V))
uh = FEFunction(V,x)

i = integrate(uh,trian,quad)

@test sum(i) + 1.0 ≈ 1.0

end # module
