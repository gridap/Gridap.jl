include("../../src/Integration/BoundaryCellQuadratures.jl")

module BoundaryCellQuadraturesTests

using Test
using Gridap
using Gridap.CellValuesGallery
using ..BoundaryCellQuadratures

partition = (2,2)
trian = CartesianTriangulation(partition)
quad = CellQuadrature(trian,order=0)

l = 3
t = HEX_AXIS
hex = Polytope((t,t,t))
cell_to_polytope = ConstantCellValue(hex,l)
facet_to_cell = CellValueFromArray([1,3,2,3])
facet_to_lfacet = CellValueFromArray([4,3,1,2])

bquad = BoundaryCellQuadrature(
  quad,cell_to_polytope,facet_to_cell,facet_to_lfacet)

@test coordinates(bquad) == coordinates(quad)

@test weights(bquad) == weights(quad)

s = coordinates_in_ref_cells(bquad)
@test s[1] == [Point(0.0,1.0,0.0),]
@test s[2] == [Point(0.0,-1.0,0.0),]
@test s[3] == [Point(0.0,0.0,-1.0),]
@test s[4] == [Point(0.0,0.0,1.0),]

end # module
