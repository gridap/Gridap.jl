module ConstrainedFESpacesTests

using Gridap

D = 2
dom = fill(4,D)
model = CartesianDiscreteModel(partition=tuple(dom...))

order = 1
diritag = "boundary"
_fespace = CLagrangianFESpace(Float64,model,order,diritag)

fixeddofs = [2,5]
fespace = ConstrainedFESpace(_fespace,fixeddofs)

ufun(x) = x[1] + x[2] + 2.0

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

bh = FEBasis(fespace)
uh = zero(fespace)

cellmat = integrate(inner(bh,bh),trian,quad)
cellvec = integrate(inner(bh,uh),trian,quad)

nfree = 7

ndiri = 16

test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

#g(x) = 2.0
#
#U = TrialFESpace(fespace,g)
#
#u(x) = x[1] + x[2] + 2.0
#
#uh = interpolate(U,u)
#
#writevtk(trian,"trian",cellfields=["uh"=>uh])




end # module
