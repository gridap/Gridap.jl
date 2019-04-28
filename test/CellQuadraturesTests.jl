module CellQuadraturesTests
##
using Test
using Numa.FieldValues
using Numa.CellValues
using Numa.Quadratures
using Numa.CellQuadratures
using Numa.Geometry
using Numa.Geometry.Cartesian

using Numa: coordinates, weights, quadrature

const D = 2

p = Point{D}(1.0,1.1)

c = [ i*p for i in 1:4 ]

w = [ i*1.0 for i in 1:4 ]

l = 10

quad = ConstantCellQuadrature(c,w,l)

q = coordinates(quad)
qw = weights(quad)

@test isa(q,CellVector{Point{D}})

@test isa(qw,CellVector{Float64})

# We can iterate coordinates and weights separately

for iq in q
  @test iq == c
end

for iw in qw
  @test iw == w
end

# or simultaneously

for (iq,iw) in zip(quad)
  @test iq == c
  @test iw == w
end
ref_quad = TensorProductQuadrature(orders=(5,4))
quad2 = ConstantCellQuadrature(ref_quad,l)

@test isa(quad2,CellQuadrature)

grid = CartesianGrid(partition=(2,3))
trian = triangulation(grid)
##

order=2
# _quadrature(celltypes(trian),order)
# function _quadrature(ct::ConstantCellValue{NTuple{Z,Int}},order) where Z
ct = celltypes(trian)
t = celldata(ct)
q = quadrature(t,order=order)
ConstantCellQuadrature(q,length(ct))

quad = quadrature(trian,order=2)
quad
typeof(t)

@test isa(quad,CellQuadrature{2})

end # module CellQuadraturesTests
