module CellFieldsBench

using Numa.FieldValues
using Numa.CellArrays
using Numa.CellFunctions

l = 1000000

include("../test/CellFunctionsTestsMocks.jl")

function doloop(x)
  for xi in x
  end
end

println("+++ CellFieldsBench ( length = $l ) +++")

tiff = inner(tbv,tbv)

print("TensorBasisBasisInner ->"); @time doloop(tiff)
print("TensorBasisBasisInner ->"); @time doloop(tiff)

vexpand = expand(vbv,sfv2)

print("VectorScalarExpand ->"); @time doloop(vexpand)
print("VectorScalarExpand ->"); @time doloop(vexpand)

include("../test/IntegrationMeshesTestsMocks.jl")
using Numa.Quadratures
using Numa.CellQuadratures

imesh = DummyIntegrationMesh2D(partition=(1000,1000))
refquad = TensorProductQuadrature(orders=(2,2))
meshcoords = cellcoordinates(imesh)
quad = ConstantCellQuadrature(refquad,length(meshcoords))
points = coordinates(quad)
phi = geomap(imesh)
basis = cellbasis(imesh)
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

kmatg = inner(valsgrad, valsgrad)

print("DummyStiffnessMatrix2D ->"); @time doloop(kmatg)
print("DummyStiffnessMatrix2D ->"); @time doloop(kmatg)

end
