module ExtrusionPolytopesTests

using Test
using Gridap.Helpers
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs

@test num_faces(VERTEX) == 1
@test num_dims(VERTEX) == 0

p = Polytope(HEX_AXIS, HEX_AXIS)
test_polytope(p,optional=true)

@test get_faces(p,0,0) == [[1], [2], [3], [4]]
@test get_faces(p,1,0) == [[1, 2], [3, 4], [1, 3], [2, 4]]
@test get_faces(p,2,0) == [[1, 2, 3, 4]]
@test get_faces(p,0,1) == [[1, 3], [1, 4], [2, 3], [2, 4]]
@test get_faces(p,1,1) == [[1], [2], [3], [4]]
@test get_faces(p,2,1) == [[1, 2, 3, 4]]
@test get_faces(p,0,2) == [[1], [1], [1], [1]]
@test get_faces(p,1,2) == [[1], [1], [1], [1]]
@test get_faces(p,2,2) == [[1]]
@test get_facedims(p) == [0,0,0,0,1,1,1,1,2]
@test get_offsets(p) == [0,4,8]
@test num_facets(p) == 4
@test num_edges(p) == 4

r = Point{2,Float64}[(1, 0), (1, 0), (0, 1), (0, 1)]
@test get_edge_tangent(p) == r

r = [
  [1, 2, 3, 4], [1, 3, 2, 4], [2, 1, 4, 3], [2, 4, 1, 3],
  [3, 1, 4, 2], [3, 4, 1, 2], [4, 2, 3, 1], [4, 3, 2, 1]]
@test Set(get_vertex_permutations(p)) == Set(r)

r = Point{2,Float64}[(0, -1), (0, 1), (-1, 0), (1, 0)]
@test get_facet_normal(p) == r

@test num_faces(p) == 9
@test num_dims(p) == 2

p = Polytope(TET_AXIS, TET_AXIS, TET_AXIS)
test_polytope(p,optional=true)

@test num_faces(p) == 15
@test num_dims(p) == 3

x = Point{3,Float64}[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
@test get_vertex_coordinates(p) == x

p = SEGMENT
test_polytope(p,optional=true)
@test get_vertex_coordinates(p) == VectorValue{1,Float64}[(0,),(1,)]
@test get_edge_tangent(p) == VectorValue{1,Float64}[(1,)]
@test get_facet_normal(p) == VectorValue{1,Float64}[(-1,),(1,)]
perm = get_vertex_permutations(p)
@test perm == [[1, 2], [2, 1]]

p = VERTEX
test_polytope(p,optional=true)
@test get_vertex_coordinates(p) == VectorValue{0,Float64}[()]
@test get_edge_tangent(p) == VectorValue{0,Float64}[]
@test get_facet_normal(p) == VectorValue{0,Float64}[]
perm = get_vertex_permutations(p)
@test perm == [[1]]

test_polytope(TRI,optional=true)
test_polytope(QUAD,optional=true)
test_polytope(TET,optional=true)
test_polytope(HEX,optional=true)

perm = get_vertex_permutations(TRI)
@test Set(perm) == Set([[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]])

perm = get_vertex_permutations(QUAD)
@test Set(perm) == Set([[1, 2, 3, 4],
                        [1, 3, 2, 4],
                        [2, 1, 4, 3],
                        [2, 4, 1, 3],
                        [3, 1, 4, 2],
                        [3, 4, 1, 2],
                        [4, 2, 3, 1],
                        [4, 3, 2, 1]])

perm = get_vertex_permutations(PYRAMID)
@test Set(perm) == Set([[1, 2, 3, 4, 5],
                        [1, 3, 2, 4, 5],
                        [2, 1, 4, 3, 5],
                        [2, 4, 1, 3, 5],
                        [3, 1, 4, 2, 5],
                        [3, 4, 1, 2, 5],
                        [4, 2, 3, 1, 5],
                        [4, 3, 2, 1, 5]])

#[1, 2, 3] ⊗ {[1, 2], [2, 1]}
#[1, 3, 2]
#[2, 1, 3]
#[2, 3, 1]
#[3, 1, 2]
#[3, 2, 1]
perm = get_vertex_permutations(WEDGE)
@test Set(perm) == Set([[1, 2, 3, 4, 5, 6],
                        [1, 3, 2, 4, 6, 5],
                        [2, 1, 3, 5, 4, 6],
                        [2, 3, 1, 5, 6, 4],
                        [3, 1, 2, 6, 4, 5],
                        [3, 2, 1, 6, 5, 4],
                        [4, 5, 6, 1, 2, 3],
                        [4, 6, 5, 1, 3, 2],
                        [5, 4, 6, 2, 1, 3],
                        [5, 6, 4, 2, 3, 1],
                        [6, 4, 5, 3, 1, 2],
                        [6, 5, 4, 3, 2, 1]])

@test num_facets(SEGMENT) == 2
@test num_facets(TRI) == 3
@test num_facets(QUAD) == 4
@test num_facets(TET) == 4
@test num_facets(HEX) == 6

@test num_vertices(PYRAMID) == 5
@test num_edges(PYRAMID) == 8
@test num_facets(PYRAMID) == 5

@test is_simplex(TRI) == true
@test is_n_cube(TRI) == false

@test is_simplex(WEDGE) == false
@test is_n_cube(WEDGE) == false

@test is_simplex(QUAD) == false
@test is_n_cube(QUAD) == true

@test is_simplex(VERTEX) == true
@test is_n_cube(VERTEX) == true

d = 1
reffaces = get_reffaces(Polytope{d},WEDGE)
iface_to_ftype = get_face_type(WEDGE,d)
@test length(reffaces) == 1
@test iface_to_ftype == [1, 1, 1, 1, 1, 1, 1, 1, 1]

expected = [[[1,2],[2,1]],[[1,2],[2,1]],[[1,2],[2,1]],[[1,2],[2,1]]]
got = get_face_vertex_permutations(QUAD,1)
@test all(Set(a) == Set(b) for (a,b) in zip(got, expected))

expected = [[[1,2,3,4],[1,3,2,4],[2,1,4,3],[2,4,1,3],[3,1,4,2],[3,4,1,2],[4,2,3,1],[4,3,2,1]]]
got = get_face_vertex_permutations(QUAD,2)
@test all(Set(a) == Set(b) for (a,b) in zip(got, expected))

expected = [
  [[1]],[[1]],[[1]],[[1]],
  [[1,2],[2,1]],[[1,2],[2,1]],[[1,2],[2,1]],[[1,2],[2,1]],
  [[1,2,3,4],[1,3,2,4],[2,1,4,3],[2,4,1,3],[3,1,4,2],[3,4,1,2],[4,2,3,1],[4,3,2,1]]]
got = get_face_vertex_permutations(QUAD)
@test all(Set(a) == Set(b) for (a,b) in zip(got, expected))

expected = [[[1,2],[2,1]],[[1,2],[2,1]],[[1,2],[2,1]]]
got = get_face_vertex_permutations(TRI,1)
@test all(Set(a) == Set(b) for (a,b) in zip(got, expected))

expected = [[[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 5, 6, 3, 4, 7, 8], [1, 3, 2, 4, 5, 7, 6, 8], [1, 3, 5, 7, 2, 4, 6, 8], [1, 5, 2, 6, 3, 7, 4, 8], [1, 5, 3, 7, 2, 6, 4, 8], [2, 1, 4, 3, 6, 5, 8, 7], [2, 1, 6, 5, 4, 3, 8, 7], [2, 4, 1, 3, 6, 8, 5, 7], [2, 4, 6, 8, 1, 3, 5, 7], [2, 6, 1, 5, 4, 8, 3, 7], [2, 6, 4, 8, 1, 5, 3, 7], [3, 1, 4, 2, 7, 5, 8, 6], [3, 1, 7, 5, 4, 2, 8, 6], [3, 4, 1, 2, 7, 8, 5, 6], [3, 4, 7, 8, 1, 2, 5, 6], [3, 7, 1, 5, 4, 8, 2, 6], [3, 7, 4, 8, 1, 5, 2, 6], [4, 2, 3, 1, 8, 6, 7, 5], [4, 2, 8, 6, 3, 1, 7, 5], [4, 3, 2, 1, 8, 7, 6, 5], [4, 3, 8, 7, 2, 1, 6, 5], [4, 8, 2, 6, 3, 7, 1, 5], [4, 8, 3, 7, 2, 6, 1, 5], [5, 1, 6, 2, 7, 3, 8, 4], [5, 1, 7, 3, 6, 2, 8, 4], [5, 6, 1, 2, 7, 8, 3, 4], [5, 6, 7, 8, 1, 2, 3, 4], [5, 7, 1, 3, 6, 8, 2, 4], [5, 7, 6, 8, 1, 3, 2, 4], [6, 2, 5, 1, 8, 4, 7, 3], [6, 2, 8, 4, 5, 1, 7, 3], [6, 5, 2, 1, 8, 7, 4, 3], [6, 5, 8, 7, 2, 1, 4, 3], [6, 8, 2, 4, 5, 7, 1, 3], [6, 8, 5, 7, 2, 4, 1, 3], [7, 3, 5, 1, 8, 4, 6, 2], [7, 3, 8, 4, 5, 1, 6, 2], [7, 5, 3, 1, 8, 6, 4, 2], [7, 5, 8, 6, 3, 1, 4, 2], [7, 8, 3, 4, 5, 6, 1, 2], [7, 8, 5, 6, 3, 4, 1, 2], [8, 4, 6, 2, 7, 3, 5, 1], [8, 4, 7, 3, 6, 2, 5, 1], [8, 6, 4, 2, 7, 5, 3, 1], [8, 6, 7, 5, 4, 2, 3, 1], [8, 7, 4, 3, 6, 5, 2, 1], [8, 7, 6, 5, 4, 3, 2, 1]]]
got = get_face_vertex_permutations(HEX,3)
@test all(Set(a) == Set(b) for (a,b) in zip(got, expected))

expected = [[[1, 2, 3, 4], [1, 2, 4, 3], [1, 3, 2, 4], [1, 3, 4, 2], [1, 4, 2, 3], [1, 4, 3, 2], [2, 1, 3, 4], [2, 1, 4, 3], [2, 3, 1, 4], [2, 3, 4, 1], [2, 4, 1, 3], [2, 4, 3, 1], [3, 1, 2, 4], [3, 1, 4, 2], [3, 2, 1, 4], [3, 2, 4, 1], [3, 4, 1, 2], [3, 4, 2, 1], [4, 1, 2, 3], [4, 1, 3, 2], [4, 2, 1, 3], [4, 2, 3, 1], [4, 3, 1, 2], [4, 3, 2, 1]]]
got = get_face_vertex_permutations(TET,3)
@test all(Set(a) == Set(b) for (a,b) in zip(got, expected))

v,p = simplexify(QUAD)
@test p == TRI

v,p = simplexify(TRI)
@test p == TRI

v,p = simplexify(HEX)
@test p == TET

v,p = simplexify(TET)
@test p == TET

v,p = simplexify(SEGMENT)
@test p == SEGMENT
@test is_simplex(p)

v,p = simplexify(VERTEX)
@test p == VERTEX
@test is_simplex(p)

for dim in 0:5
    QD = Polytope(ntuple(d->HEX_AXIS, dim))
    SD = Polytope(ntuple(d->TET_AXIS, dim))
    @test is_n_cube(QD)
    @test is_simplex(SD)

    vss,p = simplexify(QD)
    @test length(vss) == factorial(dim)
    for vs in vss
        @test allunique(vs)
        @test all(1 ≤ v ≤ 2^dim for v in vs)
    end
    @test allunique(vss)
    @test p == SD
    @test is_simplex(p)

    if dim <= 3
      vss,p = simplexify(QD,positive=true)
      @test length(vss) == factorial(dim)
      for vs in vss
          @test allunique(vs)
          @test all(1 ≤ v ≤ 2^dim for v in vs)
      end
      @test allunique(vss)
      @test p == SD
      @test is_simplex(p)
    end

    vss,p = simplexify(SD)
    @test length(vss) == 1
    for vs in vss
        @test allunique(vs)
        @test issorted(vs)
        @test all(1 ≤ v ≤ dim+1 for v in vs)
    end
    @test p == SD
    @test is_simplex(p)
end

@test SEGMENT == ExtrusionPolytope(HEX_AXIS) == ExtrusionPolytope(TET_AXIS)
@test TRI == ExtrusionPolytope(HEX_AXIS, TET_AXIS)
@test QUAD == ExtrusionPolytope(TET_AXIS, HEX_AXIS)
@test TET == ExtrusionPolytope(HEX_AXIS, TET_AXIS, TET_AXIS)
@test HEX == ExtrusionPolytope(TET_AXIS, HEX_AXIS, HEX_AXIS)
@test PYRAMID == ExtrusionPolytope(TET_AXIS, HEX_AXIS, TET_AXIS)
@test WEDGE == ExtrusionPolytope(HEX_AXIS, TET_AXIS, HEX_AXIS)

#using Gridap.Fields
#
#using .ReferenceFEs: _polytopenfaces
#using .ReferenceFEs: NFace
#using .ReferenceFEs: DFace
#using .ReferenceFEs: _polytopemesh
#using .ReferenceFEs: _dimfrom_fs_dimto_fs
#using .ReferenceFEs: _nface_to_nfacedim
#using .ReferenceFEs: _nfaces_vertices
#using .ReferenceFEs: _face_normals
#using .ReferenceFEs: _edge_tangents
#using .ReferenceFEs: _admissible_permutations
#using .ReferenceFEs: _admissible_permutations_n_cube
#
#using Profile
#using ProfileView
#
#
#extrusion = Point{0,Int}()
#p = DFace(extrusion)
#perms = _admissible_permutations(p)
#display(perms)
#
#extrusion = Point(HEX_AXIS)
#p = DFace(extrusion)
#perms = _admissible_permutations(p)
#display(perms)
#
#extrusion = Point(TET_AXIS)
#p = DFace(extrusion)
#perms = _admissible_permutations(p)
#display(perms)
#
#extrusion = Point(TET_AXIS,TET_AXIS)
#p = DFace(extrusion)
#perms = _admissible_permutations(p)
#display(perms)
#
#extrusion = Point(HEX_AXIS,HEX_AXIS)
#p = DFace(extrusion)
#perms = _admissible_permutations(p)
#display(perms)
#
#extrusion = Point(TET_AXIS, TET_AXIS, TET_AXIS)
#p = DFace(extrusion)
#perms = _admissible_permutations(p)
#display(perms)
#
#extrusion = Point(HEX_AXIS, HEX_AXIS, HEX_AXIS)
#p = DFace(extrusion)
#perms = _admissible_permutations(p)
#display(perms)
#
#
#f = DFace{2}(p,21)
#
#f = DFace{3}(p,27)


end # module
