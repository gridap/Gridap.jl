module FaceDofPermutationsTests

# The signed face DoF permutations of the moment-based reference FEs, built in
# closed form by `moment_face_own_dofs_permutations`, against a numerical
# oracle sharing none of that reasoning: an affine symmetry Φ of the *whole*
# reference cell realising the vertex permutation of a face, applied to the
# DoFs, S_{ji} = ℓ̂_j(Φ^*ψ̂_i), which must be the signed permutation matrix the
# tables describe.

using Test
using LinearAlgebra
using Gridap
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: MomentBasedDofBasis

# All affine maps x ↦ Ax + b of the reference polytope onto itself, with the
# vertex permutation each induces (vertex i goes to vertex img[i]), solved by
# least squares from the vertex coordinates for every admissible permutation.
function cell_symmetries(p::Polytope{D}) where D
  x = get_vertex_coordinates(p)
  X = hcat([v[i] for v in x, i in 1:D], ones(length(x)))
  map(get_vertex_permutations(p)) do σ
    Y = X[σ, 1:D]
    M = X \ Y
    @assert norm(X*M - Y) < 1e-12
    A = TensorValue{D,D}(Tuple(vec(Matrix(round.(Int, transpose(M[1:D, :]))))))
    b = VectorValue{D}(Tuple(round.(Int, M[D+1, :])))
    (A, b, σ)
  end
end

# ℓ̂ ∘ Φ^*, as a DoF basis: Σₙ mₙ ⊙ (Φ^*ψ)(xₙ) = Σₙ (Tᵀmₙ) ⊙ ψ(Φxₙ), with
# Φ^*ψ = T ψ∘Φ the pullback of the element's Piola map, T = Aᵀ (covariant),
# |det A| A⁻¹ (contravariant) or 1 (none).
function twist(dofs::MomentBasedDofBasis, A, b, conf)
  nodes = map(x -> A⋅x + b, dofs.nodes)
  Tᵀ = conf isa CurlConformity ? A :
       conf isa DivConformity ? abs(det(A))*transpose(inv(A)) : nothing
  f = isnothing(Tᵀ) ? identity : (m -> Tᵀ⋅m)
  moments = map(mat -> map(f, mat), dofs.face_moments)
  MomentBasedDofBasis(nodes, moments, dofs.face_nodes, dofs.face_own_moms, dofs.operator)
end

function test_face_permutations(reffe)
  p = get_polytope(reffe)
  D = num_dims(p)
  conf = Conformity(reffe)
  dofs, shapefuns = get_dof_basis(reffe), get_shapefuns(reffe)
  own = get_face_own_dofs(reffe)
  perms = get_face_own_dofs_permutations(reffe)
  signs = get_face_own_dofs_signs(reffe)
  face_vertices = reduce(vcat, [get_faces(p, d, 0) for d in 0:D])
  vperms = get_face_vertex_permutations(p)
  syms = cell_symmetries(p)
  for d in 1:D-1, f in get_dimrange(p, d)
    isempty(own[f]) && continue
    if length(perms[f]) == 1 # weights without a closed form: pindex 1 only
      @test perms[f][1] == 1:length(own[f])
      continue
    end
    @test length(perms[f]) == length(vperms[f])
    fv = face_vertices[f]
    for (pindex, π) in enumerate(vperms[f])
      A, b, _ = syms[findfirst(s -> s[3][fv] == fv[invperm(π)], syms)]
      S = evaluate(twist(dofs, A, b, conf), shapefuns)[own[f], own[f]]
      P = zeros(size(S))
      foreach(i -> P[perms[f][pindex][i], i] = signs[f][pindex][i], axes(S, 2))
      @test S ≈ P atol=1e-10
      pindex == 1 && @test P == I
    end
  end
end

# the faces (of dimension 1 to D-1) whose weights have no closed form, and so
# keep the one-entry table
function unsupported_faces(reffe)
  p = get_polytope(reffe)
  perms = get_face_own_dofs_permutations(reffe)
  [f for d in 1:num_dims(p)-1 for f in get_dimrange(p, d) if length(perms[f]) == 1]
end

for (p, order) in ((TRI, 2), (TET, 1), (QUAD, 2), (HEX, 1))
  test_face_permutations(ReferenceFE(p, raviart_thomas, order; change_dof=false))
  test_face_permutations(ReferenceFE(p, nedelec, order; change_dof=false))
end
for (p, order) in ((TRI, 2), (TET, 2))
  test_face_permutations(ReferenceFE(p, bdm, order; change_dof=false))
end
for (p, order) in ((TRI, 3), (QUAD, 3), (HEX, 2))
  test_face_permutations(ReferenceFE(p, modal_lagrangian, order; change_dof=false))
end
test_face_permutations(ReferenceFE(QUAD, sdiv, 3))
test_face_permutations(ReferenceFE(HEX, sdiv, 2))
test_face_permutations(ReferenceFE(QUAD, scurl, 3))
test_face_permutations(ReferenceFE(HEX, scurl, 3))
test_face_permutations(ReferenceFE(HEX, trimmed_sdiv, 2))
test_face_permutations(ReferenceFE(HEX, trimmed_scurl, 2))

# Nédélec's triangular face weights are 1-forms, mixed by the face symmetries:
# those faces, and only those, have no table
nd_tet = ReferenceFE(TET, nedelec, 1; change_dof=false)
@test unsupported_faces(nd_tet) == collect(get_dimrange(TET, 2))
@test isempty(unsupported_faces(ReferenceFE(HEX, nedelec, 1; change_dof=false)))
@test isempty(unsupported_faces(ReferenceFE(TET, raviart_thomas, 1; change_dof=false)))

end # module
