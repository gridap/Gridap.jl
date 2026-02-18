module GeneralPolytopesTests

using Test
using Gridap
using Gridap.ReferenceFEs

using Gridap.ReferenceFEs: check_polytope_graph
using Gridap.ReferenceFEs: simplexify_surface
using Gridap.ReferenceFEs: simplexify_interior
using Gridap.Io

p = Polygon(TRI)
test_polytope(p,optional=true)
@test check_polytope_graph(p)
@test get_faces(p,0,0) == [[1],[2],[3]]
@test get_faces(p,1,0) == [[1,2],[2,3],[3,1]]
@test get_faces(p,2,0) == [[1,2,3]]
@test get_faces(p,0,1) == [[1,3],[1,2],[2,3]]
@test get_faces(p,1,1) == [[1],[2],[3]]
@test get_faces(p,2,1) == [[1,2,3]]
@test get_faces(p,0,2) == [[1],[1],[1]]
@test get_faces(p,1,2) == [[1],[1],[1]]
@test get_faces(p,2,2) == [[1]]
@test get_facedims(p) == [0,0,0,1,1,1,2]
@test Polytope{2}(p,1) === p
@test Polytope{1}(p,1) == SEGMENT
@test Polytope{0}(p,1) == VERTEX
@test simplexify(p) == ([[1,2,3]],TRI)
@test ReferenceFEs.signed_area(p) == 0.5
quad = Quadrature(p,2)
test_quadrature(quad)

p = Polygon(QUAD)
test_polytope(p,optional=true)
@test check_polytope_graph(p)
@test get_faces(p,0,0) == [[1],[2],[3],[4]]
@test get_faces(p,1,0) == [[1,2],[2,3],[3,4],[4,1]]
@test get_faces(p,2,0) == [[1,2,3,4]]
@test get_faces(p,0,1) == [[1,4],[1,2],[2,3],[3,4]]
@test get_faces(p,1,1) == [[1],[2],[3],[4]]
@test get_faces(p,2,1) == [[1,2,3,4]]
@test get_faces(p,0,2) == [[1],[1],[1],[1]]
@test get_faces(p,1,2) == [[1],[1],[1],[1]]
@test get_faces(p,2,2) == [[1]]
@test get_facedims(p) == [0,0,0,0,1,1,1,1,2]
@test Polytope{2}(p,1) === p
@test Polytope{1}(p,1) == SEGMENT
@test Polytope{0}(p,1) == VERTEX
@test simplexify(p) == ([[1,2,3],[1,3,4]],TRI)
@test ReferenceFEs.signed_area(p) == 1.0
quad = Quadrature(p,2)
test_quadrature(quad)

p = Polyhedron(TET)
test_polytope(p,optional=true)
@test check_polytope_graph(p)
@test get_faces(p,0,0) == [[1],[2],[3],[4]]
#@test get_faces(p,1,0) == [[1,2],[1,4],[1,3],[2,3],[2,4],[3,4]]
@test get_faces(p,1,0) == [[1,3],[1,4],[1,2],[2,4],[2,3],[3,4]]
#@test get_faces(p,2,0) == [[1,2,3],[1,4,2],[1,3,4],[2,4,3]]
@test get_faces(p,2,0) == [[1,3,2], [1,4,3], [1,2,4], [2,3,4]]
@test get_faces(p,3,0) == [[1,2,3,4]]
#@test get_faces(p,0,1) == [[1,2,3],[1,4,5],[3,4,6],[2,5,6]]
@test get_faces(p,0,1) == [[1,2,3], [3,4,5], [1,5,6], [2,4,6]]
@test get_faces(p,1,1) == [[1],[2],[3],[4],[5],[6]]
#@test get_faces(p,2,1) == [[1,4,3],[2,5,1],[3,6,2],[5,6,4]]
@test get_faces(p,2,1) == [[1,5,3], [2,6,1], [3,4,2], [5,6,4]]
@test get_faces(p,3,1) == [[1,2,3,4,5,6]]
#@test get_faces(p,0,2) == [[1,2,3],[1,2,4],[1,3,4],[2,3,4]]
@test get_faces(p,0,2) == [[1, 2, 3], [1, 3, 4], [1, 2, 4], [2, 3, 4]]
#@test get_faces(p,1,2) == [[1,2],[2,3],[1,3],[1,4],[2,4],[3,4]]
@test get_faces(p,1,2) == [[1,2],[2,3],[1,3],[3,4],[1,4],[2,4]]
@test get_faces(p,2,2) == [[1],[2],[3],[4]]
@test get_faces(p,3,2) == [[1,2,3,4]]
@test get_faces(p,0,3) == [[1],[1],[1],[1]]
@test get_faces(p,1,3) == [[1],[1],[1],[1],[1],[1]]
@test get_faces(p,2,3) == [[1],[1],[1],[1]]
@test get_faces(p,3,3) == [[1]]
@test get_facedims(p) == [0,0,0,0,1,1,1,1,1,1,2,2,2,2,3]
@test Polytope{3}(p,1) === p
@test isa(Polytope{2}(p,1),Polygon)
@test Polytope{1}(p,1) == SEGMENT
@test Polytope{0}(p,1) == VERTEX
@test simplexify(p) == ([[1,2,3,4]],TET)
x,t = simplexify_surface(p)
@test t == get_faces(p,2,0)
x,t = simplexify_interior(p)
@test t == [[1,2,3,4]]
@test ReferenceFEs.signed_volume(p) == 1/6
@test ReferenceFEs.is_convex(p)
quad = Quadrature(p,2)
test_quadrature(quad)

p = Polyhedron(HEX)
test_polytope(p,optional=true)
@test check_polytope_graph(p)
@test get_faces(p,0,0) == [[1],[2],[3],[4],[5],[6],[7],[8]]
@test get_faces(p,1,0) == [
  [1,5],[1,2],[1,3],[2,6],[2,4],[3,7],[3,4],[4,8],[5,7],[5,6],[6,8],[7,8]]
@test get_faces(p,2,0) == [
  [1,5,7,3],[1,2,6,5],[1,3,4,2],[2,4,8,6],[3,7,8,4],[5,6,8,7]]
@test get_faces(p,3,0) == [[1,2,3,4,5,6,7,8]]
@test get_faces(p,0,1) == [
  [1,2,3],[2,4,5],[3,6,7],[5,7,8],[1,9,10],[4,10,11],[6,9,12],[8,11,12]]
@test get_faces(p,1,1) == [[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12]]
@test get_faces(p,2,1) == [
  [1,9,6,3],[2,4,10,1],[3,7,5,2],[5,8,11,4],[6,12,8,7],[10,11,12,9]]
@test get_faces(p,3,1) == [[1,2,3,4,5,6,7,8,9,10,11,12]]
@test get_faces(p,0,2) == [
  [1,2,3],[2,3,4],[1,3,5],[3,4,5],[1,2,6],[2,4,6],[1,5,6],[4,5,6]]
@test get_faces(p,1,2) == [
  [1,2],[2,3],[1,3],[2,4],[3,4],[1,5],[3,5],[4,5],[1,6],[2,6],[4,6],[5,6]]
@test get_faces(p,2,2) == [[1],[2],[3],[4],[5],[6]]
@test get_faces(p,3,2) == [[1,2,3,4,5,6]]
@test get_faces(p,0,3) == [[1],[1],[1],[1],[1],[1],[1],[1]]
@test get_faces(p,1,3) == [[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]]
@test get_faces(p,2,3) == [[1],[1],[1],[1],[1],[1]]
@test get_faces(p,3,3) == [[1]]
@test get_facedims(p) == [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3]
@test Polytope{3}(p,1) === p
@test isa(Polytope{2}(p,1),Polygon)
@test Polytope{1}(p,1) == SEGMENT
@test Polytope{0}(p,1) == VERTEX
@test simplexify(p) == (
  [[1,2,4,8],[1,2,8,6],[1,3,7,8],[1,3,8,4],[1,5,6,8],[1,5,8,7]],TET)
x,t = simplexify_surface(p)
@test t == [
  [1,5,7],[1,7,3],[1,2,6],[1,6,5],[1,3,4],[1,4,2],
  [2,4,8],[2,8,6],[3,7,8],[3,8,4],[5,6,8],[5,8,7]]
x,t = simplexify_interior(p)
@test t == [[1,2,4,8],[1,2,8,6],[1,3,7,8],[1,3,8,4],[1,5,6,8],[1,5,8,7]]
@test ReferenceFEs.signed_volume(p) == 1.0
@test ReferenceFEs.is_convex(p)
quad = Quadrature(p,2)
test_quadrature(quad)

p = Polygon([
  Point(0.0, 0.0), Point(3.0, 0.0), Point(3.0, 2.0), Point(2.0, 2.0), 
  Point(2.0, 1.0), Point(1.0, 1.0), Point(1.0, 2.0), Point(0.0, 2.0)
])
X = get_vertex_coordinates(p)
T, simp = simplexify(p)
polys = ReferenceFEs.convexify(p)
@test length(polys) == 3
@test all(ReferenceFEs.is_convex.(polys))
@test sum(ReferenceFEs.signed_area.(polys)) == 5.0
@test ReferenceFEs.signed_area(p) == 5.0
q = ReferenceFEs.extrude(p)

p = Polygon(map(x -> Point(x[1],x[2],0.0), get_vertex_coordinates(p)))
X = get_vertex_coordinates(p)
T, simp = simplexify(p)
polys = ReferenceFEs.convexify(p)
@test length(polys) == 3
@test all(ReferenceFEs.is_convex.(polys))
@test sum((norm∘ReferenceFEs.signed_area).(polys)) == 5.0
@test norm(ReferenceFEs.signed_area(p)) == 5.0

# IO

poly_2d = Polygon([Point(0.0,0.0), Point(1.0,0.0), Point(1.0,1.0), Point(0.0,1.0)])
dict = to_dict(poly_2d)
@test dict[:D] == 2
@test dict[:Dp] == 2
@test dict[:isopen] == false
poly_2d_r = from_dict(GeneralPolytope,dict)
@test get_vertex_coordinates(poly_2d_r) == get_vertex_coordinates(poly_2d)
@test ReferenceFEs.get_graph(poly_2d_r) == ReferenceFEs.get_graph(poly_2d)
@test isopen(poly_2d_r) == isopen(poly_2d)
@test num_dims(poly_2d_r) == num_dims(poly_2d)
@test num_point_dims(poly_2d_r) == num_point_dims(poly_2d)

tri_io = Polygon([Point(0.0,0.0), Point(1.0,0.0), Point(0.5,1.0)])
dict_tri = to_dict(tri_io)
tri_r = from_dict(GeneralPolytope,dict_tri)
@test get_vertex_coordinates(tri_r) == get_vertex_coordinates(tri_io)
@test ReferenceFEs.get_graph(tri_r) == ReferenceFEs.get_graph(tri_io)

poly_3d = Polyhedron(TET)
dict_3d = to_dict(poly_3d)
@test dict_3d[:D] == 3
@test dict_3d[:Dp] == 3
poly_3d_r = from_dict(GeneralPolytope,dict_3d)
@test get_vertex_coordinates(poly_3d_r) == get_vertex_coordinates(poly_3d)
@test ReferenceFEs.get_graph(poly_3d_r) == ReferenceFEs.get_graph(poly_3d)

poly_hex = Polyhedron(HEX)
dict_hex = to_dict(poly_hex)
poly_hex_r = from_dict(GeneralPolytope,dict_hex)
@test get_vertex_coordinates(poly_hex_r) == get_vertex_coordinates(poly_hex)
@test ReferenceFEs.get_graph(poly_hex_r) == ReferenceFEs.get_graph(poly_hex)

penta = Polygon([Point(cos(2π*i/5), sin(2π*i/5)) for i in 0:4])
dict_p = to_dict(penta)
penta_r = from_dict(GeneralPolytope,dict_p)
@test get_vertex_coordinates(penta_r) ≈ get_vertex_coordinates(penta)
@test ReferenceFEs.get_graph(penta_r) == ReferenceFEs.get_graph(penta)

json_str = to_json(poly_2d)
poly_2d_rj = from_json(GeneralPolytope,json_str)
@test get_vertex_coordinates(poly_2d_rj) == get_vertex_coordinates(poly_2d)
@test ReferenceFEs.get_graph(poly_2d_rj) == ReferenceFEs.get_graph(poly_2d)

end # module
