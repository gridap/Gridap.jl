module BubbleRefFEsTests

using Test
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.TensorValues

# Test for TRI 
p = TRI
T = Float64
reffe = BubbleRefFE(T, p)
@test reffe.polytope == p
@test reffe.ndofs == 1
face_dofs = [Int[] for _ in 1:num_faces(p)]
face_dofs[end] = [1]
@test reffe.face_dofs == face_dofs
@test typeof(reffe.prebasis) <: BubbleMonomialBasis
@test typeof(reffe.conformity) <: L2Conformity

# Test for TET (tetrahedron)
p = TET
reffe = BubbleRefFE(T, p)
@test reffe.polytope == p
@test reffe.ndofs == 1
face_dofs = [Int[] for _ in 1:num_faces(p)]
face_dofs[end] = [1]
@test reffe.face_dofs == face_dofs
@test typeof(reffe.prebasis) <: BubbleMonomialBasis
@test typeof(reffe.conformity) <: L2Conformity

# Test not implemented for unsupported polytope (e.g., QUAD)
p = QUAD
@test_throws ErrorException BubbleRefFE(T, p)

# Test for `VectorValue` type
p = TRI
D = 2
T = VectorValue{D, Float64}
reffe = BubbleRefFE(T, p)
@test reffe.polytope == p
@test reffe.ndofs == D
face_dofs = [Int[] for _ in 1:num_faces(p)]
face_dofs[end] = 1:D
@test reffe.face_dofs == face_dofs
@test typeof(reffe.prebasis) <: BubbleMonomialBasis
@test typeof(reffe.conformity) <: L2Conformity

end
