module DiscFESpacesTests

using Test
using Gridap

D = 2
dom = fill(4,D)
model = CartesianDiscreteModel(partition=tuple(dom...))

order = 1
reffe = PDiscRefFE(Float64,D,order)

fespace = DiscFESpace(reffe,model)

ufun(x) = x[1] + x[2]

trian = Triangulation(model)

nfree = 48
ndiri = 0

trian = Triangulation(model)
quad = CellQuadrature(trian,degree=2)

bh = FEBasis(fespace)
uh = zero(fespace)

cellmat = integrate(inner(bh,bh),trian,quad)
cellvec = integrate(inner(bh,uh),trian,quad)

test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

@test celldofids(fespace).data == 1:nfree

uh = interpolate(fespace,ufun)

e = ufun - uh

el2 = sqrt(sum(integrate(inner(e,e),trian,quad)))
@test el2 < 1.0e-10

#ufun(x) = sin(x[1]) * cos(x[2])
#
#uh = interpolate(fespace,ufun)
#
#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh])


end # module
