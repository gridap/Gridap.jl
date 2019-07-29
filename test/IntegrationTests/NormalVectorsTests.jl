include("../../src/Integration/NormalVectors.jl")
module NormalVectorsTests

using Test
using Gridap
using ..NormalVectors

tags = [23,24,25]
model = CartesianDiscreteModel(partition=(3,4,2))

trian = Triangulation(model)
btrian = BoundaryTriangulation(model,tags)
bquad = CellQuadrature(btrian,order=2)
q = coordinates(bquad)

phi = CellGeomap(trian)

grid = Grid(model,3)

polytopes = CellPolytopes(grid)

n = NormalVector(phi,btrian.descriptor,polytopes)

n_q = evaluate(n,q)
_ = collect(n_q)

bphi = CellGeomap(btrian)

x = evaluate(bphi,q)

writevtk(btrian,"btrian")
writevtk(x,"x",pointdata=["n"=>n_q])

end # module
