module ReferenceFEInterfacesTests

using Test
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs

# Abstract API

struct MockRefFEName <: ReferenceFEName end
@test_throws ErrorException ReferenceFE(VERTEX, MockRefFEName(), 0)

struct MockRefFE{D} <: ReferenceFE{D} end
D = 0
reffe = MockRefFE{D}()
@test_throws ErrorException num_dofs(reffe)
@test_throws ErrorException get_polytope(reffe)
@test_throws ErrorException get_prebasis(reffe)
@test_throws ErrorException get_dof_basis(reffe)
@test_throws ErrorException Conformity(reffe)
@test_throws ErrorException Conformity(reffe,nothing)
@test_throws ErrorException Conformity(reffe,:L2)
@test Conformity(reffe, L2Conformity()) == L2Conformity()
@test_throws ErrorException get_face_dofs(reffe)
@test_throws ErrorException get_face_own_dofs(reffe, nothing)

@test num_dims(reffe) == D
@test num_dims(typeof(reffe)) == D
@test num_cell_dims(reffe) == D
@test num_cell_dims(typeof(reffe)) == D
@test num_point_dims(reffe) == D
@test num_point_dims(typeof(reffe)) == D

D = 2
T = Float64
order = 1
prebasis = MonomialBasis(Val(D),T,order)

polytope = QUAD
x = get_vertex_coordinates(polytope)
dofs = LagrangianDofBasis(T,x)

ndofs = 4

face_own_dofs = [[1],[2],[3],[4],Int[],Int[],Int[],Int[],Int[]]
a = [1]
b = Int[]
face_own_dofs_permutations = [[a],[a],[a],[a],[b,b],[b,b],[b,b],[b,b],fill(b,8)]
face_dofs = [[1],[2],[3],[4],[1,2],[3,4],[1,3],[2,4],[1,2,3,4]]
metadata = nothing

struct MockConformity <: Conformity end

conf = MockConformity()
reffe = GenericRefFE{typeof(conf)}(
  ndofs, polytope, prebasis, dofs, conf, metadata ,face_dofs)

function ReferenceFEs.get_face_own_dofs(reffe::GenericRefFE{MockConformity},::MockConformity)
  face_own_dofs
end

function ReferenceFEs.get_face_own_dofs_permutations(reffe::GenericRefFE{MockConformity},::MockConformity)
  face_own_dofs_permutations
end

test_reference_fe(reffe)

@test get_face_own_dofs(reffe,L2Conformity()) == [[], [], [], [], [], [], [], [], [1, 2, 3, 4]]
@test get_face_own_dofs_permutations(reffe,L2Conformity()) == [[[]], [[]], [[]], [[]], [[]], [[]], [[]], [[]], [[1, 2, 3, 4]]]

shapefuns = get_shapefuns(reffe)

@test evaluate(shapefuns,x) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

@test get_own_dofs_permutations(reffe) == face_own_dofs_permutations[end]

# dof prebasis inversion from given shape functions
dofprebasis = dofs
shapefuns = prebasis # just to get nontrivial dofs change of basis matrix

reffe = GenericRefFE{typeof(conf)}(
  ndofs, polytope, dofprebasis, conf, metadata, face_dofs, shapefuns)

test_reference_fe(reffe)

@test get_face_own_dofs(reffe,L2Conformity()) == [[], [], [], [], [], [], [], [], [1, 2, 3, 4]]
@test get_face_own_dofs_permutations(reffe,L2Conformity()) == [[[]], [[]], [[]], [[]], [[]], [[]], [[]], [[]], [[1, 2, 3, 4]]]
@test get_own_dofs_permutations(reffe) == face_own_dofs_permutations[end]

dofs = get_dof_basis(reffe)

@test evaluate(dofs,shapefuns) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

end # module
