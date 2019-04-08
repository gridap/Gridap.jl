using Numa.Quadratures
using Numa.CellQuadratures
using Numa.CellIntegration
using Numa.CellValues
using Numa.CellFunctions
using Numa

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
fun(x::Point{2}) = x[1]*x[2] + x[1]
gradfun(x::Point{2}) = VectorValue(x[2] + 1.0, x[1])
import Numa: gradient
gradient(::typeof(fun)) = gradfun
uphys = fun ∘ phi
kmat = integrate(a(uphys,uphys),imesh,quad)
##
