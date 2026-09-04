module PΛRefFEsTests
# Reference FE layer for the rotating P_rΛ¹ and trimmed P_r⁻Λ¹ bases:
# construction, DOF/shape-function duality, identity face-own-dof
# permutations for every pindex, factories, and the geometric
# decomposition API.

using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Arrays: evaluate
using Gridap.Fields: Point
using LinearAlgebra
using Test

simplex(D) = D == 2 ? TRI : TET
make_reffe(name, D, r) = name == :rotating ? RotatingPΛRefFE(Float64, simplex(D), r) :
                                             TrimmedPΛRefFE(Float64, simplex(D), r)

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
    # shapefuns == prebasis
    pts = [Point(ntuple(_ -> rand()/(D+1), D)) for _ in 1:5]
    @test evaluate(get_shapefuns(rf), pts) == evaluate(b, pts)
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
  # factories: the T-less signature defaults to Float64
  @test ReferenceFE(TRI, rotating_pλ, 2) isa GenericRefFE{RotatingPΛName,2}
  @test ReferenceFE(TET, trimmed_pλ, 1)  isa GenericRefFE{TrimmedPΛName,3}
  @test num_dofs(ReferenceFE(TRI, rotating_pλ, 2)) ==
        num_dofs(RotatingPΛRefFE(Float64, TRI, 2))
  for (name, sing) in ((:rotating, rotating_pλ), (:trimmed, trimmed_pλ))
    rf32 = ReferenceFE(TRI, sing, Float32, 2)
    @test get_prebasis(rf32) isa PolynomialBasis{2,VectorValue{2,Float32}}
    @test num_dofs(rf32) == num_dofs(make_reffe(name, 2, 2))
    @test get_prebasis(ReferenceFE(TRI, sing, Float64, 2)) isa
          PolynomialBasis{2,VectorValue{2,Float64}}
  end
  # only defined on simplices
  @test_throws AssertionError RotatingPΛRefFE(Float64, QUAD, 1)
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
