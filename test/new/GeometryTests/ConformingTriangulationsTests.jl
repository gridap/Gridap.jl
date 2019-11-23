module ConformingTriangulationsTests

using Test
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs

import Gridap.ReferenceFEs: get_node_coordinates
import Gridap.Geometry: get_cell_nodes
import Gridap.Geometry: get_reffes
import Gridap.Geometry: get_cell_types

struct ConformingTrianMock <: ConformingTriangulation{2,2} end

function get_node_coordinates(::ConformingTrianMock)
  Point{2,Float64}[(0,0),(1,0),(2,0),(1,1),(2,1),(0,2),(2,2)]
end

function get_cell_nodes(::ConformingTrianMock)
  [[1,2,6,4],[2,3,4,5],[1,5,7],[3,7,6]]
end

function get_reffes(::ConformingTrianMock)
  order = 1
  tri3 = LagrangianRefFE(Float64,TRI,order)
  quad4 = LagrangianRefFE(Float64,QUAD,order)
  [tri3, quad4]
end

function get_cell_types(::ConformingTrianMock)
  [2,2,1,1]
end

trian = ConformingTrianMock()
test_conforming_triangulation(trian)

q1i = Point(0.25,0.25)
np1 = 3
q1 = fill(q1i,np1)
q2i = Point(0.5,0.5)
np2 = 4
q2 = fill(q2i,np2)
q = CompressedArray([q1,q2],get_cell_types(trian))

cell_map = get_cell_map(trian)
x = evaluate(cell_map,q)

x1i = Point(0.5, 0.75)
x2i = Point(1.5, 0.5)
x3i = Point(1.0, 0.75)
x4i = Point(1.5, 1.0)
x1 = fill(x1i,np2)
x2 = fill(x2i,np2)
x3 = fill(x3i,np1)
x4 = fill(x4i,np1)
x = [x1,x2,x3,x4]
test_array_of_fields(cell_map,q,x)

end # module
