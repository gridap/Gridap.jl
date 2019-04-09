##
using Numa
using Numa.Quadratures
using Numa.CellQuadratures
using Numa.CellIntegration
using Numa.CellValues
using Numa.CellFunctions
import Numa.CellIntegration: cellcoordinates, cellbasis

using Numa.CellValues: IndexCellArray
import Numa.CellValues: cellsize
using Numa.CellIntegration
import Numa.CellIntegration: cellcoordinates, cellbasis
using Numa.Polytopes
using Numa.RefFEs
using Numa.FieldValues

import Numa: gradient
include("CellIntegrationTestsMocks.jl")
##
##
polytope = Polytope(Polytopes.PointInt{2}(1,1))
reffe = LagrangianRefFE{2,ScalarValue}(polytope,[1,1])
basis = reffe.shfbasis
cellb = CellBasisFromSingleInterpolation(basis)
cellcoords = DummyCellCoordinates2D(partition=(2,2))
imesh = DummyIntegrationMesh2D(cellcoords,cellb)
refquad = TensorProductQuadrature(orders=(2,2))
meshcoords = cellcoordinates(imesh)
ncells = length(meshcoords)
# @santiagobadia : I would like a method that given a mesh, it provides the
# number of active cells
quad = ConstantCellQuadrature(refquad,ncells)
phi = geomap(imesh)
basis = cellbasis(imesh)
physbasis = attachgeomap(basis,phi)
a(v,u) = inner(∇(v),∇(u)) + inner(v,u)
V = physbasis
U = physbasis
# fun(x::Point{2}) = x[1]*x[2] + x[1]
# gradfun(x::Point{2}) = VectorValue(x[2] + 1.0, x[1])
fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
uphys = fun ∘ phi
kmat = integrate(a(uphys,uphys),imesh,quad)
sum(kmat)
kmat = integrate(a(V,U),imesh,quad)
sum(kmat)
##
