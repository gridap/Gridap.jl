module GLSRefFEsTests

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

# A *traceless* matrix field that is exactly of total degree `r`, i.e. a member
# of the element's space GLS_r(K) = P_r(K;M₀)
function poly_tensor(r)
  x -> begin
    a = 1.0 + x[1]^r
    b = 2.0 - x[2]^r
    c = 3.0 + (x[1] * x[2])^(r ÷ 2)
    TensorValue(a, b, c, -a)          # trace a + (-a) = 0
  end
end

# Generalized Vandermonde matrix Q[l,i] = ℓₗ(F*(Ψ̂ᵢ)): the physical GLS DoFs of
# the triangle `verts` evaluated on the co-contravariant Piola push-forward of
# the reference shape functions, computed from scratch by quadrature on the
# physical edges. `σ` gives, per edge, whether the cell traverses it along the
# global direction; a reversed edge flips the odd-degree Legendre weights.
# The interior DoFs are *defined* as the push-forward of the reference ones --
# they are owned by the cell and shared with nobody -- so their rows are the
# identity.
function gls_vandermonde(reffe, verts, σ; degree=12)
  p = get_polytope(reffe)
  Ψ̂ = get_shapefuns(reffe)
  ndofs = num_dofs(reffe)

  Jt = jacobian_t(verts)
  J = transpose(Jt)
  iJ = inv(J)
  detJ = det(Jt)
  F = affine_map(verts)

  own = get_face_own_dofs(reffe)
  nv = num_faces(p, 0)
  order = length(own[nv+1]) - 1
  μb = LegendreBasis(Val(1), Float64, order)

  Q = Matrix{Float64}(I, ndofs, ndofs)
  quad = Quadrature(SEGMENT, degree)
  ŝ = get_coordinates(quad)
  ŵ = get_weights(quad)

  for e in 1:num_faces(p, 1)
    v̂a, v̂b = get_face_coordinates(p, 1)[e]
    xa, xb = F(v̂a), F(v̂b)
    L = norm(xb - xa)                # physical edge length: ds = L dŝ
    t = (xb - xa) / L                # nt-moments are bilinear, both flip: either sense
    nrm = ReferenceFEs._rot90(t)

    x̂ = [v̂a + s[1] * (v̂b - v̂a) for s in ŝ]
    Ψ̂x = evaluate(Ψ̂, x̂)
    μ = evaluate(μb, ŝ)
    for (i, d) in enumerate(own[nv+e])
      parity = (σ[e] < 0 && isodd(i - 1)) ? -1.0 : 1.0
      for j in 1:ndofs
        # the mixed Piola push-forward, φ = det(J)⁻¹ J⁻ᵀ φ̂ Jᵀ
        Mtn = [(t ⋅ ((1 / detJ) * (transpose(iJ) ⋅ ψ ⋅ Jt))) ⋅ nrm for ψ in view(Ψ̂x, :, j)]
        Q[d, j] = parity * sum(ŵ .* L .* Mtn .* μ[:, i])
      end
    end
  end
  Q
end

############################################################################################
# Reference element
############################################################################################

function test_gls_reffe(r)
  reffe = GLSRefFE(Float64, TRI, r)
  ndofs = 3*(r + 1) * (r + 2) ÷ 2

  @test reffe isa GenericRefFE
  @test get_name(reffe) isa GopalakrishnanLedererSchoberl
  @test num_dofs(reffe) == ndofs
  # augmented: the prebasis is the ambient P_r(K;M), 4/3 the size of the element
  @test length(get_prebasis(reffe)) == 4 * (r + 1) * (r + 2) ÷ 2
  @test Conformity(reffe) == DivConformity()
  @test Pushforward(GopalakrishnanLedererSchoberl) == ReferenceFEs.CoContraVariantPiolaMap()
  @test reffe == ReferenceFE(TRI, gls, Float64, r)
# NB: `test_reference_fe` is not applicable here -- it asserts
# `num_dofs == length(prebasis)`, which an augmented element breaks by design.

  # r+1 DoFs per edge, the rest owned by the cell, none on the vertices
  @test length.(get_face_own_dofs(reffe)) ==
        [0, 0, 0, r + 1, r + 1, r + 1, 3*r * (r + 1) ÷ 2]

  # Unisolvency
  @test evaluate(get_dof_basis(reffe), get_shapefuns(reffe)) ≈ Matrix(I, ndofs, ndofs)
end

test_gls_reffe(0)
test_gls_reffe(1)
test_gls_reffe(2)

@test_throws ErrorException GLSRefFE(Float64, QUAD, 0)
@test_throws ErrorException GLSRefFE(Float64, TET, 0)
@test_throws ErrorException GLSRefFE(Float64, TRI, -1)

############################################################################################
# The DoFs are the claimed functionals
############################################################################################

function test_gls_dofs(r)
  reffe = GLSRefFE(Float64, TRI, r)
  prebasis = get_prebasis(reffe)
  V = evaluate(get_dof_basis(reffe), prebasis)
  own = get_face_own_dofs(reffe)
  ndofs = num_dofs(reffe)

  # The edge DoFs are ∫ₑ (t⋅Mn) μᵢ ds against the Legendre basis
  μb = LegendreBasis(Val(1), Float64, r)
  quad = Quadrature(SEGMENT, 12)
  ŝ = get_coordinates(quad)
  ŵ = get_weights(quad)
  for e in 1:3
    v̂a, v̂b = get_face_coordinates(TRI, 1)[e]
    L = norm(v̂b - v̂a)
    t = (v̂b - v̂a) / L
    n = ReferenceFEs._rot90(t)
    x̂ = [v̂a + s[1] * (v̂b - v̂a) for s in ŝ]
    φx = evaluate(prebasis, x̂)
    μ = evaluate(μb, ŝ)
    for (i, d) in enumerate(own[3+e]), j in 1:ndofs
      @test V[d, j] ≈ sum(ŵ .* L .* [(t ⋅ φ) ⋅ n for φ in view(φx, :, j)] .* μ[:, i]) atol = 1e-12
    end
  end

  # The interior DoFs are ∫_K (t_e⋅Mn_e) q dK against P_{r-1}(K), one set per edge
  if r > 0
    cb = BernsteinBasisOnSimplex(Val(2), Float64, r - 1)
    cq = Quadrature(TRI, 12)
    xc, wc = get_coordinates(cq), get_weights(cq)
    φc = evaluate(prebasis, xc)
    qc = evaluate(cb, xc)
    ts, ns = ReferenceFEs._edge_frames(TRI)
    nq = length(cb)
    for (k, d) in enumerate(own[7]), j in 1:ndofs
      e = (k - 1) ÷ nq + 1                 # which edge's frame this set uses
      i = (k - 1) % nq + 1                 # which weight
      want = sum(wc .* [(ts[e] ⋅ φ) ⋅ ns[e] for φ in view(φc, :, j)] .* qc[:, i])
      @test V[d, j] ≈ want atol = 1e-12
    end
  end
end

test_gls_dofs(0)
test_gls_dofs(1)
test_gls_dofs(2)

############################################################################################
# Transformation to a physical cell
############################################################################################

test_cells = (
  [Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0)],  # the reference cell
  [Point(0.3, 0.1), Point(1.7, 0.2), Point(0.6, 2.3)],  # det J > 0
  [Point(0.0, 0.0), Point(0.0, 1.0), Point(1.0, 0.0)],  # det J < 0
  [Point(-1.0, 0.5), Point(2.0, -0.4), Point(0.1, 0.3)],
)
test_signs = ((1.0, 1.0, 1.0), (-1.0, 1.0, 1.0), (1.0, -1.0, -1.0),
              (-1.0, -1.0, -1.0))

function test_gls_change_of_basis(r)
  reffe = GLSRefFE(Float64, TRI, r)
  ndofs = num_dofs(reffe)

  for verts in test_cells, σ in test_signs
    Jt = jacobian_t(verts)
    P = copy(evaluate(FESpaces.EdgeScalingChangeOfBasis(reffe, false), Jt, σ))
    Pinvt = copy(evaluate(FESpaces.EdgeScalingChangeOfBasis(reffe, true), Jt, σ))

    @test transpose(Pinvt) * P ≈ Matrix(I, ndofs, ndofs)

    # Unlike Morley and Argyris, this change of basis is diagonal: the
    # nt-moments are preserved up to the scalar ‖J t̂ₑ‖.
    @test P ≈ Diagonal(diag(P))

    Q = gls_vandermonde(reffe, verts, σ)
    @test Q * P ≈ Matrix(I, ndofs, ndofs)
  end

  # The edge entries really are the edge length ratio ‖J t̂ₑ‖
  verts = test_cells[2]
  Jt = jacobian_t(verts)
  P = copy(evaluate(FESpaces.EdgeScalingChangeOfBasis(reffe, false), Jt, (1.0, 1.0, 1.0)))
  own = get_face_own_dofs(reffe)
  for e in 1:3
    L = norm(get_edge_tangent(TRI)[e] ⋅ Jt)
    for d in own[3+e]
      @test P[d, d] ≈ L
    end
  end
  for d in own[7]
    @test P[d, d] ≈ 1.0
  end
end

test_gls_change_of_basis(0)
test_gls_change_of_basis(1)
test_gls_change_of_basis(2)

# At order 0 the only edge weight is the constant one, so orientation cannot
# matter at all: the change of basis ignores σ entirely.
reffe0 = GLSRefFE(Float64, TRI, 0)
Jt0 = jacobian_t(test_cells[2])
P⁺ = copy(evaluate(FESpaces.EdgeScalingChangeOfBasis(reffe0, false), Jt0, (1.0, 1.0, 1.0)))
P⁻ = copy(evaluate(FESpaces.EdgeScalingChangeOfBasis(reffe0, false), Jt0, (-1.0, -1.0, -1.0)))
@test P⁺ ≈ P⁻

############################################################################################
# FE space
#
# The discrete space must not depend on how each cell happens to list its
# vertices, so this runs on a sorted mesh and on the same mesh permuted.
############################################################################################

function test_gls_fe_space(r, model)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 2*r + 6)
  V = FESpace(model, ReferenceFE(gls, Float64, r))

  topo = get_grid_topology(model)
  @test num_free_dofs(V) == (r + 1) * num_faces(topo, 1) +
                            (3*r * (r + 1) ÷ 2) * num_cells(model)

  # A global traceless P_r matrix field is interpolated exactly
  τ = poly_tensor(r)
  τh = interpolate(τ, V)
  @test sqrt(sum(∫((τh - τ) ⊙ (τh - τ))dΩ)) < 1e-12

  # A field outside the space, to test conformity
  w(x) = TensorValue(sin(2*x[1]), cos(3*x[2]), sin(x[1] + x[2]), -sin(2*x[1]))
  wh = interpolate(w, V)
  @test sqrt(sum(∫((wh - w) ⊙ (wh - w))dΩ)) > 1e-6

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 2*r + 6)
  n = get_normal_vector(Λ).⁺
  t = Operation(ReferenceFEs._rot90)(n)   # the unit tangent of the shared edge

  # Only the normal-tangential component is continuous -- the defining property
  @test sqrt(sum(∫((t ⋅ jump(wh) ⋅ n) * (t ⋅ jump(wh) ⋅ n))dΛ)) < 1e-12
  @test sqrt(sum(∫(jump(wh) ⊙ jump(wh))dΛ)) > 1e-3
end

sorted3 = unit_square(3)
permuted3 = permute_cells(unit_square(3))

test_gls_fe_space(0, sorted3)
test_gls_fe_space(1, sorted3)
test_gls_fe_space(2, sorted3)
test_gls_fe_space(0, permuted3)
test_gls_fe_space(1, permuted3)
test_gls_fe_space(2, permuted3)

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
