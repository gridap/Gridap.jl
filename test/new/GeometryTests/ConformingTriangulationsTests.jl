module ConformingTriangulationsTests

using Test
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs

using Gridap.Geometry: ConformingTrianMock

trian = ConformingTrianMock()
test_conforming_triangulation(trian)

q1i = Point(0.5,0.5)
np1 = 4
q1 = fill(q1i,np1)

q2i = Point(0.25,0.25)
np2 = 3
q2 = fill(q2i,np2)

q = CompressedArray([q1,q2],get_cell_types(trian))

cell_map = get_cell_map(trian)
x = evaluate(cell_map,q)

x1i = Point(1.0, 0.5)
x2i = Point(1.5, 0.25)
x3i = Point(1.5, 0.75)
x4i = Point(1.25, 1.0)
x5i = Point(0.5, 0.75)
x1 = fill(x1i,np1)
x2 = fill(x2i,np2)
x3 = fill(x3i,np2)
x4 = fill(x4i,np1)
x5 = fill(x5i,np1)
x = [x1,x2,x3,x4,x5]

test_array_of_fields(cell_map,q,x)

end # module
