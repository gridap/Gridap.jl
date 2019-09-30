module CellQuadraturesTests

using Test
using Gridap

p = Point(1.0,1.1)
l = 10

ref_quad = TensorProductQuadrature(orders=(5,4))
quad2 = ConstantCellQuadrature(ref_quad,l)

coo_cell = coordinates(quad2)
wei_cell = weights(quad2)

coo = coordinates(ref_quad)
wei = weights(ref_quad)

for i in 1:quad2.length
  @test coo_cell[1] == coo
  @test wei_cell[1] == wei
end

@test isa(quad2,CellQuadrature)

@test string(quad2) == "CellQuadrature object"
s = "CellQuadrature object:\n dim: 2\n ncells: 10"
@test sprint(show,"text/plain",quad2) == s

end # module CellQuadraturesTests
