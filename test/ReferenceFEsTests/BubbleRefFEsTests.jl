module BubbleRefFEsTests

using Test
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.TensorValues
using Gridap.Fields: LinearCombinationFieldVector
using Statistics: mean
using Gridap.Arrays: evaluate

# mini bubble tests
et = Float64
for p in [SEGMENT, TRI, QUAD, TET, HEX]
	for T in [et, VectorValue{2, et}, VectorValue{3, et}]
		reffe = BubbleRefFE(T, p)
		test_reference_fe(reffe)

		N = num_components(T)
		@test num_dofs(reffe) == N
		@test Conformity(reffe) == L2Conformity()
		@test get_polytope(reffe) == p

		face_dofs = fill(Int[], num_faces(p))
		face_dofs[end] = 1:N
		@test get_face_dofs(reffe) == face_dofs

		prebasis = get_prebasis(reffe)
		@test length(prebasis) == N
		@test prebasis isa LinearCombinationFieldVector

		shapefuns = get_shapefuns(reffe)
		@test length(shapefuns) == N

		dofs = get_dof_basis(reffe)
		xs = get_face_coordinates(p)
		bxs = map(mean, xs)
		bx0 = bxs[end]
		# dof is at the barycenter of the polytope
		for dof in dofs
			@test bx0 == dof.point
		end
		# equal to 1 at the barycenter
		val = evaluate(dofs, shapefuns)
		@test val == one(val)

		# equal to 0 at the barycenters of each d < D faces
		vals = evaluate(shapefuns, bxs[1:(end-1)])
		@test all(vals .== zero(T))
	end
end

# specific tests for TRI
p = TRI
reffe = BubbleRefFE(Float64, p)
shapefuns = get_shapefuns(reffe)
xs = [VectorValue(x, y) for x in 0:0.1:1, y in 0:0.1:1 if x + y <= 1.0]
foo((x, y)) = 27*x*y*(1-x-y)

_is_approx(x, y) = isapprox(x, y, atol = 1e-14)
@test all(map(_is_approx, evaluate(shapefuns, xs), foo.(xs)))

# specific tests for QUAD
p = QUAD
reffe = BubbleRefFE(Float64, p)
shapefuns = get_shapefuns(reffe)
xs = vec([VectorValue(x, y) for x in 0:0.1:1, y in 0:0.1:1])
foo((x, y)) = 16*x*(1-x)*y*(1-y)

@test all(map(_is_approx, evaluate(shapefuns, xs), foo.(xs)))

end
