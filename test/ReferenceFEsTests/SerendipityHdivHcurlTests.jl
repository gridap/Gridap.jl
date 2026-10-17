module SerendipityHdivHcurlTests

using Test
using LinearAlgebra
using Gridap.Helpers
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs

# [Arnold & Awanou, Math. Comp. 83 (2014), (22)]
aa_dim(D,r,k) = sum(2^(D-d)*binomial(D,d)*binomial(r-d+2k,d)*binomial(d,k) for d in k:min(D,r÷2+k))
# [Gillette & Kloefkorn, Math. Comp. 88 (2019), Table 1], order k = 0, 1, …
const GK_DIMS = Dict((2,1) => [4,10,17,26], (3,1) => [12,36,66,111], (3,2) => [6,21,45,82])

const ELEMENTS = (
  # name, constructor, conformity, dimension, orders on QUAD and on HEX
  (sdiv,          SerendipityDivRefFE,         DivConformity(),  (D,o) -> aa_dim(D,o,D-1),       1:4, 1:3),
  (scurl,         SerendipityCurlRefFE,        CurlConformity(), (D,o) -> aa_dim(D,o,1),         1:4, 1:3),
  (trimmed_sdiv,  TrimmedSerendipityDivRefFE,  DivConformity(),  (D,o) -> GK_DIMS[(D,D-1)][o+1], 0:3, 0:2),
  (trimmed_scurl, TrimmedSerendipityCurlRefFE, CurlConformity(), (D,o) -> GK_DIMS[(D,1)][o+1],   0:3, 0:3),
)

for (name, RefFE, conf, dim, qorders, horders) in ELEMENTS, (p, orders) in ((QUAD,qorders), (HEX,horders)), o in orders
  D = num_dims(p)
  reffe = RefFE(Float64, p, o)
  @test reffe isa GenericRefFE{typeof(name)}
  @test num_dofs(reffe) == dim(D,o)
  @test length(get_prebasis(reffe)) == num_dofs(reffe)
  @test Conformity(reffe) == conf
  @test all(isempty, get_face_own_dofs(reffe)[get_dimrange(p,0)])

  dofs = get_dof_basis(reffe)
  @test cond(evaluate(dofs, get_prebasis(reffe))) < 1e3
  @test evaluate(dofs, get_shapefuns(reffe)) ≈ I atol=1e-12
end

@test ReferenceFE(QUAD, sdiv, Float64, 2) isa GenericRefFE{SerendipityDiv}
@test ReferenceFE(HEX, scurl, Float64, 1) isa GenericRefFE{SerendipityCurl}
@test ReferenceFE(QUAD, trimmed_sdiv, Float64, 0) isa GenericRefFE{TrimmedSerendipityDiv}
@test ReferenceFE(HEX, trimmed_scurl, Float64, 1) isa GenericRefFE{TrimmedSerendipityCurl}
@test Pushforward(SerendipityDiv) == ContraVariantPiolaMap()
@test Pushforward(TrimmedSerendipityCurl) == CoVariantPiolaMap()
for (_, RefFE) in ELEMENTS
  @test_throws ErrorException RefFE(Float64, TRI, 1)
  @test_throws ErrorException RefFE(Float64, QUAD, -1)
end

# Ownership: facet moments against 𝓟ᵣ(f), tangential face moments in 3D
for r in 1:3
  reffe = SerendipityDivRefFE(Float64, HEX, r)
  own = get_face_own_dofs(reffe)
  @test all(isempty, own[get_dimrange(HEX,1)])
  @test all(f -> length(f) == binomial(r+2,2), own[get_dimrange(HEX,2)])
  @test length(own[end]) == 3*binomial(r+1,3)
  reffe = SerendipityCurlRefFE(Float64, HEX, r)
  own = get_face_own_dofs(reffe)
  @test all(e -> length(e) == r+1, own[get_dimrange(HEX,1)])
  @test all(f -> length(f) == 2*binomial(r,2), own[get_dimrange(HEX,2)])
end

# The lowest orders are Raviart-Thomas and Nédélec
for p in (QUAD,HEX)
  x = [Point(ntuple(i -> (0.1+0.2*i+0.13*j)%1, num_dims(p))) for j in 1:10]
  span(b) = (v = evaluate(b,x); reduce(hcat, [[c for w in col for c in Tuple(w)] for col in eachcol(v)]))
  same(A,B) = rank([A B]; rtol=1e-10) == rank(A; rtol=1e-10) == rank(B; rtol=1e-10)
  @test same(span(get_shapefuns(TrimmedSerendipityDivRefFE(Float64,p,0))),
             span(get_shapefuns(ReferenceFE(p,raviart_thomas,Float64,0))))
  @test same(span(get_shapefuns(TrimmedSerendipityCurlRefFE(Float64,p,0))),
             span(get_shapefuns(ReferenceFE(p,nedelec,Float64,0))))
end

# The interior DoFs of the trimmed H(curl) element pair with curl 𝓗ₖ, not ∇𝓗ₖ
# [Gillette & Kloefkorn, (4.1)]: with gradients, unisolvence is lost at every
# even k.
for k in 1:4
  reffe = TrimmedSerendipityCurlRefFE(Float64, QUAD, k)
  pre = get_prebasis(reffe)
  emom(φ,μ,ds) = Broadcasting(Operation(*))(Broadcasting(Operation(⋅))(φ,get_edge_tangent(ds)),μ)
  cmom(φ,μ,ds) = Broadcasting(Operation(⋅))(φ,μ)
  moments = Tuple[(get_dimrange(QUAD,1), emom, LegendreBasis(Val(1),Float64,k,Polynomials._p_filter))]
  if k > 1
    push!(moments, (get_dimrange(QUAD,2), cmom, LegendreBasis(Val(2),VectorValue{2,Float64},k-2,Polynomials._p_filter)))
  end
  gradients = Polynomials._exact_poly_basis(Legendre,Float64,2,k-1,1)
  V = evaluate(MomentBasedDofBasis(QUAD, pre, [moments; (get_dimrange(QUAD,2), cmom, gradients)]), pre)
  @test (rank(V) == size(V,1)) == isodd(k)
end

# The face DoFs are signed permutations of the global ones: checked against the
# reference DoFs taken through each symmetry Φ(x) = Ax+b of the cube realising a
# face permutation, ℓ̂ⱼ(Φ*ψ̂ᵢ), with Gridap's Piola maps.
function twist(dofs::MomentBasedDofBasis, A, b)
  nodes = map(x -> A⋅x + b, dofs.nodes)
  moments = map(m -> map(mᵢ -> A⋅mᵢ, m), dofs.face_moments)
  MomentBasedDofBasis(nodes, moments, dofs.face_nodes, dofs.face_own_moms, dofs.operator)
end

function test_face_permutations(reffe)
  p = get_polytope(reffe)
  D = num_dims(p)
  dofs, shapefuns = get_dof_basis(reffe), get_shapefuns(reffe)
  own = get_face_own_dofs(reffe)
  perms = ReferenceFEs.get_metadata(reffe).face_own_dofs_permutations
  signs = ReferenceFEs.get_metadata(reffe).face_own_dofs_signs
  @test get_face_own_dofs_permutations(reffe, Conformity(reffe)) === perms
  @test get_face_own_dofs_permutations(reffe, L2Conformity()) != perms

  face_vertices = reduce(vcat, [get_faces(p,d,0) for d in 0:D])
  vperms = get_face_vertex_permutations(p)
  syms = ReferenceFEs._cube_symmetries(p)
  for d in 1:D-1, f in get_dimrange(p,d)
    isempty(own[f]) && continue
    fv = face_vertices[f]
    for (pindex, π) in enumerate(vperms[f])
      A, b, _ = syms[findfirst(s -> s[3][fv] == fv[invperm(π)], syms)]
      S = evaluate(twist(dofs,A,b), shapefuns)[own[f],own[f]]
      P = zeros(size(S))
      foreach(i -> P[perms[f][pindex][i],i] = signs[f][pindex][i], axes(S,2))
      @test S ≈ P atol=1e-10
    end
  end
end

test_face_permutations(SerendipityDivRefFE(Float64,QUAD,3))
test_face_permutations(SerendipityDivRefFE(Float64,HEX,2))
test_face_permutations(SerendipityCurlRefFE(Float64,QUAD,3))
test_face_permutations(SerendipityCurlRefFE(Float64,HEX,3))
test_face_permutations(TrimmedSerendipityDivRefFE(Float64,HEX,2))
test_face_permutations(TrimmedSerendipityCurlRefFE(Float64,HEX,2))

end # module
