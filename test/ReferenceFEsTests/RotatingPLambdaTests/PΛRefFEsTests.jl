module PΛRefFEsTests
# Reference FE layer for the rotating P_rΛ¹ and trimmed P_r⁻Λ¹ bases:
# construction, DOF/shape-function duality, identity face-own-dof
# permutations for every pindex, factories, and the geometric
# decomposition API.

using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Arrays: evaluate
using LinearAlgebra
using Test

make_reffe(name, D, r) = name == :rotating ? RotatingPΛRefFE(D, r) : TrimmedPΛRefFE(D, r)

# ─── Construction, duality, identity dof perms ───────────────────────────────

@testset "PΛ reffe construction, duality, identity dof perms" begin
  for name in (:rotating, :trimmed), (D, r) in ((2,1), (2,2), (2,3), (3,1), (3,2))
    rf = make_reffe(name, D, r)
    @test rf isa GenericRefFE
    @test get_name(rf) isa (name == :rotating ? RotatingPΛName : TrimmedPΛName)
    @test Conformity(rf) == CurlConformity()
    b = get_prebasis(rf)
    @test num_dofs(rf) == length(b)
    # exact duality: the dof basis is the dual basis of the shape functions
    E = evaluate(get_dof_basis(rf), get_shapefuns(rf))
    @test maximum(abs, E - I) < 1e-12
    # identity face-own-dof permutations for EVERY pindex of every face
    poly  = get_polytope(rf)
    perms = get_face_own_dofs_permutations(rf, Conformity(rf))
    own   = get_face_own_dofs(rf)
    vtx_perms = ReferenceFEs.get_face_vertex_permutations(poly)
    @test length(perms) == num_faces(poly)
    for gf in 1:num_faces(poly)
      @test length(perms[gf]) == length(vtx_perms[gf])
      @test all(p == collect(1:length(own[gf])) for p in perms[gf])
    end
  end
  # factories
  @test ReferenceFE(TRI, rotating_pλ, 2) isa GenericRefFE{RotatingPΛName,2}
  @test ReferenceFE(TET, trimmed_pλ, 1)  isa GenericRefFE{TrimmedPΛName,3}
  @test num_dofs(ReferenceFE(TRI, rotating_pλ, 2)) == num_dofs(RotatingPΛRefFE(2, 2))
end

# ─── Geometric decomposition API ─────────────────────────────────────────────

@testset "has_geometric_decomposition / get_face_own_funs" begin
  for name in (:rotating, :trimmed), (D, r) in ((2,2), (3,1))
    rf = make_reffe(name, D, r)
    b  = get_prebasis(rf)
    p  = get_polytope(rf)
    @test ReferenceFEs.has_geometric_decomposition(b, p, CurlConformity())
    @test ReferenceFEs.has_geometric_decomposition(b, p, L2Conformity())
    @test !ReferenceFEs.has_geometric_decomposition(b, p, ReferenceFEs.GradConformity())
    @test ReferenceFEs.get_face_own_funs(b, p, CurlConformity()) == get_face_own_dofs(rf)
  end
end

end # module
