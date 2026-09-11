module ArnoldWintherCRefFEsTests

using LinearAlgebra
using Test

using Gridap
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData

############################################################################################
# Helpers
############################################################################################

# Transposed Jacobian of the affine map from the reference triangle onto the
# triangle with vertices `verts`, in Gridap's convention Jt[i,j] = ∂Fⱼ/∂ξⁱ.
function jacobian_t(verts)
  e1 = verts[2] - verts[1]
  e2 = verts[3] - verts[1]
  J = TensorValue(e1[1], e1[2], e2[1], e2[2])  # columns e1, e2
  transpose(J)
end

affine_map(verts) = ξ -> verts[1] + (verts[2] - verts[1]) * ξ[1] +
                         (verts[3] - verts[1]) * ξ[2]

# A rebuilt model carries no face labels, so re-create the usual ones.
function tag_boundary!(model)
  topo = get_grid_topology(model)
  labels = get_face_labeling(model)
  D = num_cell_dims(model)
  for d in 0:D-1
    get_face_entity(labels, d) .= ifelse.(Geometry.get_isboundary_face(topo, d), 2, 1)
  end
  get_face_entity(labels, D) .= 1
  add_tag!(labels, "interior", [1])
  add_tag!(labels, "boundary", [2])
  model
end

# Rebuild `model` with the vertices of each cell listed in a permuted order,
# cycling through `perms`. The mesh is geometrically identical; only the local
# vertex numbering -- and hence the direction in which each cell traverses its
# edges -- changes.
function permute_cells(model, perms=([1, 2, 3], [2, 3, 1], [3, 1, 2],
                                     [2, 1, 3], [1, 3, 2], [3, 2, 1]))
  grid = get_grid(model)
  cell_nodes = get_cell_node_ids(grid)
  permuted = [collect(cell_nodes[i])[perms[mod1(i, length(perms))]]
              for i in 1:length(cell_nodes)]
  pgrid = UnstructuredGrid(get_node_coordinates(grid), Table(permuted),
                           get_reffes(grid), get_cell_type(grid))
  tag_boundary!(UnstructuredDiscreteModel(pgrid))
end

unit_square(n) = simplexify(CartesianDiscreteModel((0, 1, 0, 1), (n, n)))

# A basis of P₂(K) ∩ P₁(K)^⊥ by explicit L² projection. `ArnoldWintherCRefFEs.jl`
# takes the degree-2 members of a `DubinerBasis` instead, so this is the test's
# own oracle, sharing no code with the element it checks.
function awc_div_complement(p::Polytope{2}, x, w)
  p1 = evaluate(MonomialBasis(Val(2), Float64, 1, Polynomials._p_filter), x)
  p2 = evaluate(MonomialBasis(Val(2), Float64, 2, Polynomials._p_filter), x)
  @assert size(p1, 2) == 3 && size(p2, 2) == 6

  G = transpose(p1) * (w .* p1)                          # Gram matrix of P₁
  R = p2 - p1 * (G \ (transpose(p1) * (w .* p2)))        # the P₁-orthogonal parts

  # keep an independent set: the range of R is exactly P₂ ∩ P₁^⊥, of dimension 3
  M = Symmetric(transpose(R) * (w .* R))
  F = eigen(M)
  keep = findall(λ -> λ > 1e-10*maximum(F.values), F.values)
  @assert length(keep) == 3
  R * F.vectors[:, keep]
end

# Generalized Vandermonde matrix Q[l,i] = ℓₗ(F*(Ψ̂ᵢ)): the physical AWc DoFs of
# the triangle `verts` -- vertex components and the edge moments -- evaluated on
# the double contravariant Piola push-forward of the reference shape functions.
# The interior DoFs are *defined* as the push-forward of the reference ones, so
# their rows are the identity.
function awc_vandermonde(reffe, verts, σ; degree=12)
  p = get_polytope(reffe)
  Ψ̂ = get_shapefuns(reffe)
  ndofs = num_dofs(reffe)
  Jt = jacobian_t(verts)
  J = transpose(Jt)
  detJ = det(Jt)
  F = affine_map(verts)
  own = get_face_own_dofs(reffe)
  nv = num_faces(p, 0)
  push_(ψ) = (1 / detJ^2) * (J ⋅ ψ ⋅ transpose(J))

  Q = Matrix{Float64}(I, ndofs, ndofs)
  Ψ̂v = evaluate(Ψ̂, get_vertex_coordinates(p))
  for v in 1:nv, j in 1:ndofs
    τ = push_(Ψ̂v[v, j])
    Q[own[v][1], j] = τ[1, 1]
    Q[own[v][2], j] = τ[1, 2]
    Q[own[v][3], j] = τ[2, 2]
  end

  μb = LegendreBasis(Val(1), Float64, 1)
  quad = Quadrature(SEGMENT, degree)
  ŝ = get_coordinates(quad)
  ŵ = get_weights(quad)
  for e in 1:num_faces(p, 1)
    v̂a, v̂b = get_face_coordinates(p, 1)[e]
    xa, xb = F(v̂a), F(v̂b)
    L = norm(xb - xa)
    t = (xb - xa) / L
    n = ReferenceFEs._rot90(t)
    x̂ = [v̂a + s[1] * (v̂b - v̂a) for s in ŝ]
    Ψ̂x = evaluate(Ψ̂, x̂)
    μ = evaluate(μb, ŝ)
    dofs = own[nv+e]
    nm = length(dofs) ÷ 2
    for i in 1:nm
      par = (σ[e] < 0 && isodd(i - 1)) ? -1.0 : 1.0
      for j in 1:ndofs
        Φ = [push_(ψ) for ψ in view(Ψ̂x, :, j)]
        Q[dofs[i], j] = par * sum(ŵ .* L .* [(n ⋅ φ) ⋅ n for φ in Φ] .* μ[:, i])
        Q[dofs[nm+i], j] = par * sum(ŵ .* L .* [(n ⋅ φ) ⋅ t for φ in Φ] .* μ[:, i])
      end
    end
  end
  Q
end

############################################################################################
# Reference element
############################################################################################

reffe = ArnoldWintherCRefFE(Float64, TRI)

@test reffe isa GenericRefFE
@test get_name(reffe) isa ArnoldWintherC
@test num_dofs(reffe) == 24
# the prebasis is the ambient P₃(K;S), not the element's own 24-dimensional
# space: the augmented construction restricts to the DoFs through the shape
# functions, so `num_dofs != length(prebasis)` here by design
@test length(get_prebasis(reffe)) == 30
@test Conformity(reffe) == DivConformity()
@test Pushforward(ArnoldWintherC) == ReferenceFEs.DoubleContraVariantPiolaMap()
@test reffe == ReferenceFE(TRI, aw_c, Float64)
# NB: `test_reference_fe` is not applicable here -- it asserts
# `num_dofs == length(prebasis)`, which an augmented element breaks by design.

# three DoFs per vertex, four per edge, three interior
@test length.(get_face_own_dofs(reffe)) == [3, 3, 3, 4, 4, 4, 3]

@test evaluate(get_dof_basis(reffe), get_shapefuns(reffe)) ≈ Matrix(I, 24, 24)

@test_throws ErrorException ArnoldWintherCRefFE(Float64, QUAD)
@test_throws ErrorException ArnoldWintherCRefFE(Float64, TET)
@test_throws ErrorException ReferenceFE(TRI, aw_c, Float64, 2)

############################################################################################
# The space is what the constraints say it is
############################################################################################

Ψ = get_shapefuns(reffe)
quad_c = Quadrature(TRI, 10)
xq, wq = get_coordinates(quad_c), get_weights(quad_c)

# div Ψ ∈ P₁(K;R²): its moments against a basis of P₂ ∩ P₁^⊥ vanish
dv = evaluate(Broadcasting(divergence)(Ψ), xq)
qc = awc_div_complement(TRI, xq, wq)
for j in 1:24, c in 1:2, m in 1:3
  @test abs(sum(wq .* [d[c] for d in view(dv, :, j)] .* qc[:, m])) < 1e-12
end

# P₂(K;S) ⊂ AWc(K): the element reproduces any quadratic symmetric tensor. This
# is the check that catches a wrong constraint set -- `div Ψ ∈ P₁` and dim 24
# alone did not.
Ψxq = evaluate(Ψ, xq)
for τ in (x -> SymTensorValue{2,Float64}(1.0 + x[1]^2, 3.0 - x[1] * x[2], 2.0 + x[2]^2),
          x -> SymTensorValue{2,Float64}(x[1] * x[2], x[2]^2, x[1]^2),
          x -> SymTensorValue{2,Float64}(1.0, x[1], x[2]))
  d = evaluate(get_dof_basis(reffe), GenericField(τ))
  rec = [sum(d[j] * Ψxq[i, j] for j in 1:24) for i in eachindex(xq)]
  @test maximum(norm(rec[i] - τ(xq[i])) for i in eachindex(xq)) < 1e-11
end

############################################################################################
# Transformation to a physical cell
############################################################################################

test_cells = (
  [Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0)],
  [Point(0.3, 0.1), Point(1.7, 0.2), Point(0.6, 2.3)],
  [Point(0.0, 0.0), Point(0.0, 1.0), Point(1.0, 0.0)],
  [Point(-1.0, 0.5), Point(2.0, -0.4), Point(0.1, 0.3)],
)
test_signs = ((1.0, 1.0, 1.0), (-1.0, 1.0, 1.0), (1.0, -1.0, -1.0),
              (-1.0, -1.0, -1.0))

for verts in test_cells, σ in test_signs
  Jt = jacobian_t(verts)
  P = copy(evaluate(FESpaces.AWCChangeOfBasis(reffe, false), Jt, σ))
  Pinvt = copy(evaluate(FESpaces.AWCChangeOfBasis(reffe, true), Jt, σ))

  @test transpose(Pinvt) * P ≈ Matrix(I, 24, 24)
  @test isapprox(awc_vandermonde(reffe, verts, σ) * P, Matrix(I, 24, 24); atol=1e-10)

  if verts in (test_cells[2], test_cells[4])
    @test !isapprox(awc_vandermonde(reffe, verts, σ), Matrix(I, 24, 24))
  end
end

# The vertex block is det(J)⁻² times the matrix of H ↦ J H Jᵀ
for verts in test_cells
  Jt = jacobian_t(verts)
  J = transpose(Jt)
  B = FESpaces._congruence_matrix(J) / det(Jt)^2
  for (col, Ĥ) in enumerate((TensorValue(1.0, 0.0, 0.0, 0.0),
                             TensorValue(0.0, 1.0, 1.0, 0.0),
                             TensorValue(0.0, 0.0, 0.0, 1.0)))
    H = (J ⋅ Ĥ ⋅ transpose(J)) / det(Jt)^2
    e = VectorValue(ntuple(i -> i == col ? 1.0 : 0.0, 3))
    @test B ⋅ e ≈ VectorValue(H[1, 1], H[1, 2], H[2, 2])
  end
end

############################################################################################
# FE space
############################################################################################

function test_awc_fe_space(model)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 12)
  V = FESpace(model, ReferenceFE(TRI, aw_c, Float64))

  topo = get_grid_topology(model)
  @test num_free_dofs(V) == 3*num_faces(topo, 0) + 4*num_faces(topo, 1) +
                            3*num_cells(model)

  # P₂(K;S) ⊂ AWc, so a global quadratic symmetric tensor is interpolated exactly
  τ(x) = SymTensorValue{2,Float64}(1.0 + x[1]^2, 3.0 - x[1] * x[2], 2.0 + x[2]^2 - x[1])
  τh = interpolate(τ, V)
  @test sqrt(sum(∫((τh - τ) ⊙ (τh - τ))dΩ)) < 1e-11

  w(x) = SymTensorValue{2,Float64}(sin(2*x[1]), cos(3*x[2]), sin(x[1] + x[2]))
  wh = interpolate(w, V)
  @test sqrt(sum(∫((wh - w) ⊙ (wh - w))dΩ)) > 1e-8

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 12)
  n = get_normal_vector(Λ).⁺

  # AWc is conforming in H(div;S): the whole traction is continuous ...
  @test sqrt(sum(∫((jump(wh) ⋅ n) ⋅ (jump(wh) ⋅ n))dΛ)) < 1e-11
  # ... while the tensor itself is not
  @test sqrt(sum(∫(jump(wh) ⊙ jump(wh))dΛ)) > 1e-4
end

test_awc_fe_space(unit_square(3))
test_awc_fe_space(permute_cells(unit_square(3)))

############################################################################################
# The DIV operator
############################################################################################

div_model = unit_square(3)
Ωd = Triangulation(div_model)
dΩd = Measure(Ωd, 10)
dωd = Measure(Ωd, 10, integration_domain_style=ReferenceDomain())

Vσ = FESpace(div_model, ReferenceFE(TRI, aw_c, Float64))
Vu = FESpace(div_model, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1); conformity=:L2)

# a stress field inside the space, whose divergence is known in closed form
τd(x) = SymTensorValue{2,Float64}(1.0 + x[1]^2, 3.0 - x[1] * x[2], 2.0 + x[2]^2 - x[1])
divτd(x) = VectorValue(x[1], x[2])
σh = interpolate(τd, Vσ)
@test sqrt(sum(∫((σh - τd) ⊙ (σh - τd))dΩd)) < 1e-11    # τ really is in the space
vh = interpolate(x -> VectorValue(sin(2*x[1]), cos(3*x[2])), Vu)

@test sum(∫( DIV(σh) ⋅ vh )dωd) ≈ sum(∫( divτd ⋅ vh )dΩd)
@test sum(∫( DIV(σh) ⋅ DIV(σh) )dωd) > 1e-3             # not trivially zero

############################################################################################
# Edge orientation
############################################################################################

sorted_model = unit_square(3)
topo = get_grid_topology(sorted_model)
@test all(issorted, Geometry.get_faces(topo, 2, 0))
@test all(σ -> all(isequal(1.0), σ), FESpaces._edge_signs(sorted_model, TRI))

pmodel = permute_cells(sorted_model)
cell_vertices = Geometry.get_faces(get_grid_topology(pmodel), 2, 0)
@test !all(issorted, cell_vertices)
signs = FESpaces._edge_signs(pmodel, TRI)
@test any(σ -> any(isequal(-1.0), σ), signs)

end # module
