module CellFieldsBench

using Gridap
using Gridap.FieldValues
using Gridap.CellValues
using Gridap.CellFunctions
using Gridap.Geometry
using Gridap.Geometry.Cartesian

l = 1000000

include("../test/CellFunctionsTestsMocks.jl")

function doloop(x)
  for xi in x
  end
end

println("+++ CellFieldsBench ( length = $l ) +++")

tiff = varinner(tbv,tbv)

print("TensorBasisBasisInner ->"); @time doloop(tiff)
print("TensorBasisBasisInner ->"); @time doloop(tiff)

vexpand = expand(vbv,sfv2)

print("VectorScalarExpand ->"); @time doloop(vexpand)
print("VectorScalarExpand ->"); @time doloop(vexpand)

using Gridap.Geometry
using Gridap.Quadratures
using Gridap.CellQuadratures
using Gridap.CellIntegration

grid = CartesianGrid(partition=(1000,1000))
trian = Triangulation(grid)

meshcoords = cellcoordinates(trian)
quad = CellQuadrature(trian,order=2)
points = coordinates(quad)
phi = geomap(trian)
basis = CellBasis(trian)
physbasis = attachgeomap(basis,phi)
physbasisgrad = gradient(physbasis)
valsgrad = evaluate(physbasisgrad,points)
xg = evaluate(phi,points)
vals = evaluate(basis,points)

print("DummyMeshCoords2D ->"); @time doloop(meshcoords)
print("DummyMeshCoords2D ->"); @time doloop(meshcoords)

print("DummyCellBasis2D ->"); @time doloop(vals)
print("DummyCellBasis2D ->"); @time doloop(vals)

print("xg ->"); @time doloop(xg)
print("xg ->"); @time doloop(xg)

print("PhysBasisGradVals ->"); @time doloop(valsgrad)
print("PhysBasisGradVals ->"); @time doloop(valsgrad)

kmatg = varinner(valsgrad, valsgrad)

print("DummyStiffnessMatrix2DAtQPoints ->"); @time doloop(kmatg)
print("DummyStiffnessMatrix2DAtQPoints ->"); @time doloop(kmatg)

a(v,u) = inner(∇(v),∇(u)) + inner(v,u)
V = physbasis
U = physbasis
kmat = integrate(a(V,U),trian,quad)

print("DummyStiffnessMatrix2D ->"); @time doloop(kmat)
print("DummyStiffnessMatrix2D ->"); @time doloop(kmat)

end
