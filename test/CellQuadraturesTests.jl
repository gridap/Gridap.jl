module CellQuadraturesTests
##
using Test
using Gridap.FieldValues
using Gridap.CellValues
using Gridap.CellValues.ConstantCellValues
using Gridap.Quadratures
using Gridap.CellQuadratures
using Gridap.Geometry
using Gridap.Geometry.Cartesian

using Gridap: coordinates, weights
D = 2

p = Point{D}(1.0,1.1)
l = 10
##

# c = [ i*p for i in 1:4 ]
#
# w = [ i*1.0 for i in 1:4 ]
#
#
# quad = ConstantCellQuadrature(c,w,l)
#
# q = coordinates(quad)
# qw = weights(quad)
#
# @test isa(q,CellVector{Point{D}})
#
# @test isa(qw,CellVector{Float64})
##
# We can iterate coordinates and weights separately
#
# for iq in q
#   @test iq == c
# end
#
# for iw in qw
#   @test iw == w
# end
#
# or simultaneously

# for (iq,iw) in zip(quad)
#   @test iq == c
#   @test iw == w
# end

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

grid = CartesianGrid(partition=(2,3))
trian = Triangulation(grid)
##

order=2
# _quadrature(celltypes(trian),order)
# function _quadrature(ct::ConstantCellValue{NTuple{Z,Int}},order) where Z
ct = celltypes(trian)
t = celldata(ct)
q = Quadrature(t,order=order)
ConstantCellQuadrature(q,length(ct))

quad = CellQuadrature(trian,order=2)
quad
typeof(t)

@test isa(quad,CellQuadrature{2})

end # module CellQuadraturesTests
