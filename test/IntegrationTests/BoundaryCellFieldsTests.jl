module BoundaryCellFieldsTests

using Test
using Gridap

tags = [23,24,25]
model = CartesianDiscreteModel(partition=(3,4,2))

trian = Triangulation(model)
btrian = BoundaryTriangulation(model,tags)

ufun(x) = 2*x[1] + x[2]

u = CellField(trian,ufun)

bu = restrict(u,btrian)

quad = CellQuadrature(btrian,order=8)
q = coordinates(quad)

bphi = CellGeomap(btrian)

x = evaluate(bphi,q)

buq = evaluate(bu,q)

phi = CellGeomap(trian)
bphi2 = restrict(phi,btrian)

x2 = evaluate(bphi2,q)

@test x ≈ x2

#writevtk(btrian,"btrian",cellfields=["u"=>bu])

#writevtk(x,"x",pointdata=["u"=>buq])

end # module
