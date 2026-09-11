module ArgyrisRefFEsTests

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

# Generalized Vandermonde matrix Q[l,i] = ℓₗ(F*(Ψ̂ᵢ)) of [Kirby, (3.10)]: the
# physical Argyris DoFs of the triangle `verts` -- vertex values, gradients and
# Hessians in the global Cartesian frame, and edge normal-derivative moments --
# evaluated on the pullback of the reference shape functions. Computed from
# scratch, independently of the change-of-basis code.
function argyris_vandermonde(reffe, verts; degree=8)
  p = get_polytope(reffe)
  Ψ̂ = get_shapefuns(reffe)
  ndofs = num_dofs(reffe)

  Jt = jacobian_t(verts)
  K = inv(Jt)                 # J⁻ᵀ: ∇u = K∇û and D²u = K D²û Kᵀ
  F = affine_map(verts)
  V̂ = get_vertex_coordinates(p)

  Q = zeros(Float64, ndofs, ndofs)
  Q[1:3, :] .= evaluate(Ψ̂, V̂)

  ∇Ψ̂v = evaluate(Broadcasting(∇)(Ψ̂), V̂)
  HΨ̂v = evaluate(Broadcasting(∇∇)(Ψ̂), V̂)
  for v in 1:3, j in 1:ndofs
    g = K ⋅ ∇Ψ̂v[v, j]
    Q[3+2*v-1, j] = g[1]
    Q[3+2*v, j] = g[2]
    H = K ⋅ HΨ̂v[v, j] ⋅ transpose(K)
    Q[9+3*v-2, j] = H[1, 1]
    Q[9+3*v-1, j] = H[1, 2]
    Q[9+3*v, j] = H[2, 2]
  end

  quad = Quadrature(SEGMENT, degree)
  ŝ = get_coordinates(quad)
  ŵ = get_weights(quad)
  for k in 1:num_faces(p, 1)
    v̂a, v̂b = get_face_coordinates(p, 1)[k]
    xa, xb = F(v̂a), F(v̂b)
    L = norm(xb - xa)
    n = ReferenceFEs._rot90((xb - xa) / L)

    x̂ = [v̂a + si[1] * (v̂b - v̂a) for si in ŝ]
    ∇Ψ̂ = evaluate(Broadcasting(∇)(Ψ̂), x̂)
    for j in 1:ndofs
      Q[18+k, j] = sum(ŵ .* L .* [(K ⋅ g) ⋅ n for g in view(∇Ψ̂, :, j)])
    end
  end
  Q
end

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
# edges -- changes. Both orientations occur, since `perms` holds even and odd
# permutations.
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

############################################################################################
# Reference element
############################################################################################

reffe = ArgyrisRefFE(Float64, TRI)

@test reffe isa GenericRefFE
@test get_name(reffe) isa Argyris
@test num_dofs(reffe) == 21
@test length(get_prebasis(reffe)) == 21
@test Conformity(reffe) == H1Conformity()
@test Pushforward(Argyris) == IdentityPiolaMap()
@test reffe == ReferenceFE(TRI, argyris, Float64)
test_reference_fe(reffe)

# Six DoFs per vertex, one per edge; the cell owns none
@test get_face_own_dofs(reffe) == [[1, 4, 5, 10, 11, 12], [2, 6, 7, 13, 14, 15],
                                   [3, 8, 9, 16, 17, 18], [19], [20], [21], Int[]]

# Unisolvency: the moment form of the edge DoF gives a valid element
E = evaluate(get_dof_basis(reffe), get_shapefuns(reffe))
@test E ≈ Matrix(I, 21, 21)
@test cond(evaluate(get_dof_basis(reffe), get_prebasis(reffe))) < 1e5

# The element is defined on triangles only, at order 5 only
@test_throws ErrorException ArgyrisRefFE(Float64, QUAD)
@test_throws ErrorException ArgyrisRefFE(Float64, TET)
@test_throws ErrorException ReferenceFE(TRI, argyris, Float64, 2)

############################################################################################
# The DoFs are the claimed functionals
############################################################################################

prebasis = get_prebasis(reffe)
V = evaluate(get_dof_basis(reffe), prebasis)
V̂ = get_vertex_coordinates(TRI)

# Vertex values, then gradient components, then Hessian components -- all point
# evaluations, expressed as moments over the 0-dimensional faces
@test V[1:3, :] ≈ evaluate(prebasis, V̂)
∇φ = evaluate(Broadcasting(∇)(prebasis), V̂)
Hφ = evaluate(Broadcasting(∇∇)(prebasis), V̂)
for v in 1:3
  @test V[3+2*v-1, :] ≈ [g[1] for g in view(∇φ, v, :)]
  @test V[3+2*v, :] ≈ [g[2] for g in view(∇φ, v, :)]
  @test V[9+3*v-2, :] ≈ [h[1, 1] for h in view(Hφ, v, :)]
  @test V[9+3*v-1, :] ≈ [h[1, 2] for h in view(Hφ, v, :)]
  @test V[9+3*v, :] ≈ [h[2, 2] for h in view(Hφ, v, :)]
end

quad = Quadrature(SEGMENT, 8)
ŝ = get_coordinates(quad)
ŵ = get_weights(quad)

# Edge normal-derivative moments ∫ₑ ∇u⋅n ds. atol: most of these are exactly
# zero, and `≈` alone compares two different roundings of zero.
for k in 1:3
  v̂a, v̂b = get_face_coordinates(TRI, 1)[k]
  L = norm(v̂b - v̂a)
  n = ReferenceFEs._rot90((v̂b - v̂a) / L)
  x̂ = [v̂a + si[1] * (v̂b - v̂a) for si in ŝ]
  ∇e = evaluate(Broadcasting(∇)(prebasis), x̂)
  for j in 1:21
    @test V[18+k, j] ≈ sum(ŵ .* L .* [g ⋅ n for g in view(∇e, :, j)]) atol = 1e-12
  end
end

# The identity that collapses Kirby's node completion for the moment form of the
# edge DoF: ∫ₑ ∇u⋅t ds = u(v_b) - u(v_a), exactly, for any smooth u
Ψ = get_shapefuns(reffe)
for k in 1:3
  v̂a, v̂b = get_face_coordinates(TRI, 1)[k]
  L = norm(v̂b - v̂a)
  t = (v̂b - v̂a) / L
  x̂ = [v̂a + si[1] * (v̂b - v̂a) for si in ŝ]
  ∇Ψ = evaluate(Broadcasting(∇)(Ψ), x̂)
  Ψa, Ψb = evaluate(Ψ, [v̂a]), evaluate(Ψ, [v̂b])
  for j in 1:21
    @test sum(ŵ .* L .* [g ⋅ t for g in view(∇Ψ, :, j)]) ≈ Ψb[1, j] - Ψa[1, j] atol = 1e-12
  end
end

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

for verts in test_cells, σ in test_signs
  Jt = jacobian_t(verts)
  P = copy(evaluate(FESpaces.ArgyrisChangeOfBasis(TRI, false), Jt, σ))
  Pinvt = copy(evaluate(FESpaces.ArgyrisChangeOfBasis(TRI, true), Jt, σ))

  @test transpose(Pinvt) * P ≈ Matrix(I, 21, 21)

  Dσ = Diagonal([ones(18); collect(σ)])
  Q = Dσ * argyris_vandermonde(reffe, verts)
  @test Q * P ≈ Matrix(I, 21, 21)

  if verts in (test_cells[2], test_cells[4])
    @test !isapprox(Q, Matrix(I, 21, 21))
  end
end

# The Hessian block is the matrix of H ↦ K H Kᵀ on the components (11,12,22)
for verts in test_cells
  K = inv(jacobian_t(verts))
  B = FESpaces._congruence_matrix(K)
  for (col, Ĥ) in enumerate((TensorValue(1.0, 0.0, 0.0, 0.0),
                             TensorValue(0.0, 1.0, 1.0, 0.0),
                             TensorValue(0.0, 0.0, 0.0, 1.0)))
    H = K ⋅ Ĥ ⋅ transpose(K)
    e = VectorValue(ntuple(i -> i == col ? 1.0 : 0.0, 3))
    @test B ⋅ e ≈ VectorValue(H[1, 1], H[1, 2], H[2, 2])
  end
end

############################################################################################
# FE space
#
# The discrete space must not depend on how each cell happens to list its
# vertices, so this runs on a mesh whose cells are sorted by global vertex id and
# on the same mesh with the cell vertices permuted.
############################################################################################

function test_argyris_fe_space(model)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 12)
  V = FESpace(model, ReferenceFE(argyris, Float64, 5))

  topo = get_grid_topology(model)
  @test num_free_dofs(V) == 6*num_faces(topo, 0) + num_faces(topo, 1)

  # The space is P₅ on each cell and the DoFs are unisolvent, so a global quintic
  # must be interpolated exactly.
  u(x) = 1 + x[1] - 2*x[2] + x[1]^2 * x[2]^3 - 0.5*x[1]^5 + x[1] * x[2]^4 + 3*x[1]^3 * x[2]
  uh = interpolate(u, V)
  @test sqrt(sum(∫((uh - u) * (uh - u))dΩ)) < 1e-13

  # A field outside the space, to test conformity
  w(x) = sin(2*x[1]) * cos(3*x[2])
  wh = interpolate(w, V)
  @test sqrt(sum(∫((wh - w) * (wh - w))dΩ)) > 1e-8

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 12)

  # Argyris is C¹: both the value and the full gradient are continuous ...
  @test sqrt(sum(∫(jump(wh) * jump(wh))dΛ)) < 1e-12
  @test sqrt(sum(∫(jump(∇(wh)) ⋅ jump(∇(wh)))dΛ)) < 1e-11
  # ... but not C²
  @test sqrt(sum(∫(jump(∇∇(wh)) ⊙ jump(∇∇(wh)))dΛ)) > 1e-4
end

test_argyris_fe_space(unit_square(3))
test_argyris_fe_space(permute_cells(unit_square(3)))

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

# A sign is -1 exactly where the cell lists the edge against the order the mesh
# stores it in. The master is the mesh's own order, not `sortperm` of the global
# ids: they agree on a sorted mesh but not after permuting, differing by a global
# sign per edge.
signs = FESpaces._edge_signs(pmodel, TRI)
edge_lvertices = get_faces(TRI, 1, 0)
ptopo = get_grid_topology(pmodel)
edge_vertices = Geometry.get_faces(ptopo, 1, 0)
cell_edges = Geometry.get_faces(ptopo, 2, 1)
for cell in 1:length(cell_vertices)
  v = cell_vertices[cell]
  for (e, lv) in enumerate(edge_lvertices)
    aligned = collect(v[lv]) == collect(edge_vertices[cell_edges[cell][e]])
    @test signs[cell][e] == (aligned ? 1.0 : -1.0)
  end
end
@test any(σ -> any(isequal(-1.0), σ), signs)

end # module
