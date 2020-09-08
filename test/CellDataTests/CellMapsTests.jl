module CellMapsTests

using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Integration
using Gridap.CellData
using Gridap.Fields
using Gridap.TensorValues
using LinearAlgebra
using FillArrays
using Test

node_to_x = 2*Point{2,Float64}[
  Point(0,0),
  Point(1,0),
  Point(2,0),
  Point(0,1),
  Point(1,1),
  Point(2,1),
  Point(0,2),
  Point(1,2),
  Point(2,2)]

cell_to_node =[
  [1,2,4,5],
  [2,3,5,6],
  [4,5,7,8],
  [5,6,8,9]]

cell_to_x = [ node_to_x[nodes] for nodes in cell_to_node]

l = length(cell_to_node)

reffe = QUAD4

degree = 2
refqua = Quadrature(get_polytope(reffe),degree)
q = GenericCellPoint(Fill(get_coordinates(refqua),l))
w = Fill(get_weights(refqua),l)
npoin = num_points(refqua)

test_cell_point(q)

# cell map
ndofs = num_dofs(reffe)
s = Fill(get_shapefuns(reffe),l)
ϕ = GenericCellMap(lincomb(s,cell_to_x))
test_cell_map(ϕ)

x = ϕ(q)
test_cell_point(x)

ϕinv = inverse_map(ϕ)

end # module
