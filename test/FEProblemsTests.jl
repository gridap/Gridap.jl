using Numa.Quadratures
using Numa.CellQuadratures
using Numa.CellIntegration
using Numa.CellValues
using Numa.CellFunctions
using Numa.Commons

include("CellIntegrationTestsMocks.jl")

##
imesh = DummyIntegrationMesh2D(partition=(2,2))
refquad = TensorProductQuadrature(orders=(2,2))
meshcoords = cellcoordinates(imesh)
quad = ConstantCellQuadrature(refquad,length(meshcoords))
points = coordinates(quad)
phi = geomap(imesh)
basis = cellbasis(imesh)
physbasis = attachgeomap(basis,phi)
a(v,u) = inner(∇(v),∇(u)) + inner(v,u)
V = physbasis
U = physbasis
kmat = integrate(a(V,U),imesh,quad)
##
