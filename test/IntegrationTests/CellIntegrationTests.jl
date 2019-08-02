module CellIntegrationTests

using Test
using Gridap

partition = (3,3)
trian = CartesianTriangulation(partition)

basis = CellBasis(trian)

quad = CellQuadrature(trian,order=2)

m = varinner(basis,basis)

mmat = integrate(m,trian,quad)

@test isa(mmat,CellArray{Float64,2})

ufun(x) = 1.0

cellvol = integrate(ufun,trian,quad)

@test isa(cellvol,CellValue{Float64})

for vi in cellvol
  @assert vi ≈ (2.0/3)^2
end


# Test for simplices

model = CartesianDiscreteModel(domain=(0,2,0,2),partition=(2,2))
model = simplexify(model)
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)
cellvol = integrate(ufun,trian,quad)
vol = sum(cellvol)
@test vol ≈ 4

model = CartesianDiscreteModel(domain=(0,2,0,2,0,2),partition=(2,2,3))
model = simplexify(model)
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)
cellvol = integrate(ufun,trian,quad)
vol = sum(cellvol)
@test vol ≈ 8

end # module
