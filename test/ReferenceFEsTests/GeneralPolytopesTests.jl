module GeneralPolytopesTests

using Test
using Gridap
using Gridap.ReferenceFEs

using Gridap.ReferenceFEs: check_polytope_graph
using Gridap.ReferenceFEs: simplexify_surface
using Gridap.ReferenceFEs: simplexify_interior

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

p = Polyhedron(TET)
test_polytope(p,optional=true)
@test check_polytope_graph(p)
@test get_faces(p,0,0) == [[1],[2],[3],[4]]
@test get_faces(p,1,0) == [[1,2],[1,4],[1,3],[2,3],[2,4],[3,4]]
@test get_faces(p,2,0) == [[1,2,3],[1,4,2],[1,3,4],[2,4,3]]
@test get_faces(p,3,0) == [[1,2,3,4]]
@test get_faces(p,0,1) == [[1,2,3],[1,4,5],[3,4,6],[2,5,6]]
@test get_faces(p,1,1) == [[1],[2],[3],[4],[5],[6]]
@test get_faces(p,2,1) == [[1,4,3],[2,5,1],[3,6,2],[5,6,4]]
@test get_faces(p,3,1) == [[1,2,3,4,5,6]]
@test get_faces(p,0,2) == [[1,2,3],[1,2,4],[1,3,4],[2,3,4]]
@test get_faces(p,1,2) == [[1,2],[2,3],[1,3],[1,4],[2,4],[3,4]]
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
@test simplexify(p) == ([[1,2,4,3]],TET)
x,t = simplexify_surface(p)
@test t == get_faces(p,2,0)
x,t = simplexify_interior(p)
@test t == [[1,2,4,3]]

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
@test t == [
  [1,2,4,8],[1,2,8,6],[1,3,7,8],[1,3,8,4],[1,5,6,8],[1,5,8,7]]

end # module
