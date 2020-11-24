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

end # module
