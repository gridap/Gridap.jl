module NormalVectorsTests

using Test
using Gridap


tags = [23,24,25]
model = CartesianDiscreteModel(partition=(3,4,2))

btrian = BoundaryTriangulation(model,tags)

n = NormalVector(btrian)

sr = "NormalVector object:\n length: 20"
@test sprint(show,"text/plain",n) == sr

@test string(n) == "NormalVector object"

bquad = CellQuadrature(btrian,degree=2)
q = coordinates(bquad)

n_q = evaluate(n,q)
a = collect(n_q)

@test a[1] ≈ VectorValue{3,Float64}[(3.1225e-17, -1.0, -0.0), (0.0, -1.0, -1.11022e-16), (4.16334e-17, -1.0, -0.0), (-0.0, -1.0, -0.0)]
@test a[14] ≈ VectorValue{3,Float64}[(-1.04083e-17, -1.0, -1.38778e-17), (-0.0, -1.0, -0.0), (-0.0, -1.0, -0.0), (-0.0, -1.0, -0.0)]

bphi = CellGeomap(btrian)

x = evaluate(bphi,q)

#writevtk(btrian,"btrian")
#writevtk(x,"x",pointdata=["n"=>n_q])

strian = SkeletonTriangulation(model,"interior")
squad = CellQuadrature(strian,degree=2)
qs = coordinates(squad)

sphi = CellGeomap(strian)

n = NormalVector(strian)
x = evaluate(sphi,qs)
n_q = evaluate(n,qs)

#writevtk(strian,"strian")
#writevtk(x,"x",pointdata=["n"=>n_q])



end # module
