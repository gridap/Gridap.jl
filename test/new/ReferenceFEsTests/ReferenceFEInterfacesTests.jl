module ReferenceFEInterfacesTests

using Test
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs

D = 2
T = Float64
order = 1
prebasis = MonomialBasis{D}(T,order)

polytope = QUAD
x = get_vertex_coordinates(polytope)
dofs = LagrangianDofBasis(T,x)

facedofids = [[1],[2],[3],[4],Int[],Int[],Int[],Int[],Int[]]

reffe = GenericRefFE(polytope,prebasis,dofs,facedofids)
test_reference_fe(reffe)

shapefuns = get_shapefuns(reffe)

@test evaluate(shapefuns,x) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

reffe = LagrangianRefFE(Float64,TRI,3)

face_to_dofs = get_face_dofs(reffe)

@test face_to_dofs == [
  [1], [2], [3], [1, 2, 4, 5], [1, 3, 6, 7], [2, 3, 8, 9], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]

face_to_nodes = get_face_nodes(reffe)
@test face_to_nodes == face_to_dofs

@test has_straight_faces(reffe) == false

@test has_straight_faces(TRI3) == true

@test is_affine(TRI3) == true

@test is_affine(QUAD4) == false

end # module
