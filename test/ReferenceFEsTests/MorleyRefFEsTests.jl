module MorleyRefFEsTests

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
# physical Morley DoFs of the triangle `verts` -- vertex values and edge
# normal-derivative moments -- evaluated on the pullback of the reference shape
# functions. Computed from scratch, independently of the change-of-basis code.
function physical_vandermonde(reffe, verts; degree=6)
  p = get_polytope(reffe)
  Ψ̂ = get_shapefuns(reffe)
  ndofs = num_dofs(reffe)

  iJt = inv(jacobian_t(verts))   # ∇u = J⁻ᵀ∇û, i.e. `iJt ⋅ ∇û`
  F = affine_map(verts)

  Q = zeros(Float64, ndofs, ndofs)
  Q[1:3, :] .= evaluate(Ψ̂, get_vertex_coordinates(p))

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
      Q[3+k, j] = sum(ŵ .* L .* [(iJt ⋅ g) ⋅ n for g in view(∇Ψ̂, :, j)])
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

reffe = MorleyRefFE(Float64, TRI)

@test reffe isa GenericRefFE
@test get_name(reffe) isa Morley
@test num_dofs(reffe) == 6
@test length(get_prebasis(reffe)) == 6
@test Conformity(reffe) == H1Conformity()
@test Pushforward(Morley) == IdentityPiolaMap()
@test reffe == ReferenceFE(TRI, morley, Float64)
test_reference_fe(reffe)

# One DoF per vertex and one per edge; the cell owns none
@test get_face_own_dofs(reffe) == [[1], [2], [3], [4], [5], [6], Int[]]

# Duality of the DoF and shape function bases
E = evaluate(get_dof_basis(reffe), get_shapefuns(reffe))
@test E ≈ Matrix(I, 6, 6)

# The element is defined on triangles only, at order 2 only
@test_throws ErrorException MorleyRefFE(Float64, QUAD)
@test_throws ErrorException MorleyRefFE(Float64, TET)
@test_throws ErrorException ReferenceFE(TRI, morley, Float64, 1)

############################################################################################
# The DoFs are the claimed functionals
############################################################################################

Ψ = get_shapefuns(reffe)
prebasis = get_prebasis(reffe)
V = evaluate(get_dof_basis(reffe), prebasis)

# The first three DoFs are the vertex values ...
@test V[1:3, :] ≈ evaluate(prebasis, get_vertex_coordinates(TRI))

quad = Quadrature(SEGMENT, 6)
ŝ = get_coordinates(quad)
ŵ = get_weights(quad)

# ... and the last three the edge normal-derivative moments ∫ₑ ∇u⋅n ds. atol:
# most of these are exactly zero, and `≈` alone compares two roundings of zero.
for k in 1:3
  v̂a, v̂b = get_face_coordinates(TRI, 1)[k]
  L = norm(v̂b - v̂a)
  n = ReferenceFEs._rot90((v̂b - v̂a) / L)
  x̂ = [v̂a + si[1] * (v̂b - v̂a) for si in ŝ]
  ∇φ = evaluate(Broadcasting(∇)(prebasis), x̂)
  for j in 1:6
    @test V[3+k, j] ≈ sum(ŵ .* L .* [g ⋅ n for g in view(∇φ, :, j)]) atol = 1e-12
  end
end

# The identity that collapses Kirby's node completion for this element:
# ∫ₑ ∇u⋅t ds = u(v_b) - u(v_a), exactly, for any smooth u.
for k in 1:3
  v̂a, v̂b = get_face_coordinates(TRI, 1)[k]
  L = norm(v̂b - v̂a)
  t = (v̂b - v̂a) / L
  x̂ = [v̂a + si[1] * (v̂b - v̂a) for si in ŝ]
  ∇Ψ = evaluate(Broadcasting(∇)(Ψ), x̂)
  Ψa = evaluate(Ψ, [v̂a])
  Ψb = evaluate(Ψ, [v̂b])
  for j in 1:6
    moment = sum(ŵ .* L .* [g ⋅ t for g in view(∇Ψ, :, j)])
    @test moment ≈ Ψb[1, j] - Ψa[1, j] atol = 1e-13
  end
end

############################################################################################
# Transformation to a physical cell
############################################################################################

ts, ns = ReferenceFEs._edge_frames(TRI)

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
  P = copy(evaluate(FESpaces.MorleyChangeOfBasis(TRI, false), Jt, σ))
  Pinvt = copy(evaluate(FESpaces.MorleyChangeOfBasis(TRI, true), Jt, σ))

  # The closed forms are consistent: Pinvt is the transposed inverse of P
  @test transpose(Pinvt) * P ≈ Matrix(I, 6, 6)

  # The physical nodal basis is Ψ = linear_combination(P, F*(Ψ̂)), i.e. the
  # physical DoFs are dual to it: Q P = I with Q the Vandermonde matrix.
  # Reversing the global direction of an edge flips its DoF.
  Dσ = Diagonal([1.0, 1.0, 1.0, σ...])
  Q = Dσ * physical_vandermonde(reffe, verts)
  @test Q * P ≈ Matrix(I, 6, 6)

  # Without the change of basis, the pulled-back basis is NOT nodal. Only
  # asserted for the two generic triangles: on the reference cell, and on its
  # mirror image with every edge flipped, Q is the identity on the nose.
  if verts in (test_cells[2], test_cells[4])
    @test !isapprox(Q, Matrix(I, 6, 6))
  end
end

# The closed forms Aₖ = det(J) n̂ᵀGn̂, Bₖ = det(J) t̂ᵀGn̂ agree with a direct
# evaluation of ‖Jt̂‖ (n̂ or t̂)⋅(J⁻¹n), which is how they were derived
for verts in test_cells
  Jt = jacobian_t(verts)
  J = transpose(Jt)
  G = inv(Jt ⋅ transpose(Jt))
  for e in 1:3
    t̂, n̂ = ts[e], ns[e]
    Jt̂ = t̂ ⋅ Jt                       # = J t̂
    Jinv_n = inv(J) ⋅ ReferenceFEs._rot90(Jt̂ / norm(Jt̂))
    @test det(Jt) * (n̂ ⋅ (G ⋅ n̂)) ≈ norm(Jt̂) * (n̂ ⋅ Jinv_n)
    @test det(Jt) * (t̂ ⋅ (G ⋅ n̂)) ≈ norm(Jt̂) * (t̂ ⋅ Jinv_n)
  end
end

############################################################################################
# FE space
#
# The discrete space must not depend on how each cell happens to list its
# vertices, so this runs on a mesh whose cells are sorted by global vertex id and
# on the same mesh with the cell vertices permuted.
############################################################################################

function test_morley_fe_space(model)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 6)
  V = FESpace(model, ReferenceFE(morley, Float64, 2))

  topo = get_grid_topology(model)
  @test num_free_dofs(V) == num_faces(topo, 0) + num_faces(topo, 1)

  # The space is P₂ on each cell and the DoFs are unisolvent, so a global
  # quadratic must be interpolated exactly.
  u(x) = 1.0 + 2*x[1] - x[2] + 3*x[1]^2 - x[1] * x[2] + 0.5*x[2]^2
  uh = interpolate(u, V)
  @test sqrt(sum(∫((uh - u) * (uh - u))dΩ)) < 1e-13

  # A field outside the space, to test conformity
  w(x) = sin(2*x[1]) * cos(3*x[2])
  wh = interpolate(w, V)
  @test sqrt(sum(∫((wh - w) * (wh - w))dΩ)) > 1e-6

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 6)
  n_Λ = get_normal_vector(Λ)

  # Morley is neither C⁰ ...
  @test sqrt(sum(∫(jump(wh) * jump(wh))dΛ)) > 1e-6
  # ... nor C¹; only the facet mean of the normal derivative is continuous,
  # which is exactly the DoF the two cells share
  @test sqrt(sum(∫(jump(∇(wh) ⋅ n_Λ) * jump(∇(wh) ⋅ n_Λ))dΛ)) > 1e-6
  @test maximum(abs, get_array(∫(jump(∇(wh) ⋅ n_Λ))dΛ)) < 1e-13

  # NB: the interpolant of a non-polynomial function is NOT identical on the two
  # meshes -- see the discretized moments below. What is invariant is the space,
  # and hence the solution of a problem posed in it.
end

test_morley_fe_space(unit_square(4))
test_morley_fe_space(permute_cells(unit_square(4)))

############################################################################################
# Discretized moments
#
# `MomentBasedDofBasis` builds *discretized* moments: the integral is replaced by
# a quadrature exact on the prebasis (here of degree order + 1 = 3). The edge DoF
# is therefore the exact ∫ₑ ∇u⋅n ds only for u polynomial enough, and merely an
# approximation of it otherwise. Unlike a Piola-mapped element -- where the
# reference and physical moments are related by an exact change of variables, so
# the same quadrature computes both -- Morley's change of basis mixes the edge
# and vertex DoFs with cell-dependent coefficients, so that quadrature error is
# cell dependent. Consequence: interpolating a non-polynomial function gives
# slightly different results depending on how the cells list their vertices. The
# FE space itself, and any solution computed in it, are unaffected.
############################################################################################

dofs = get_dof_basis(reffe)

fine_quad = Quadrature(SEGMENT, 20)
ŝ_f = get_coordinates(fine_quad)
ŵ_f = get_weights(fine_quad)
function exact_edge_dofs(∇u)
  [begin
     v̂a, v̂b = get_face_coordinates(TRI, 1)[k]
     L = norm(v̂b - v̂a)
     n = ReferenceFEs._rot90((v̂b - v̂a) / L)
     sum(ŵ_f .* L .* [∇u(v̂a + si[1] * (v̂b - v̂a)) ⋅ n for si in ŝ_f])
   end for k in 1:3]
end

# exact for a quadratic ...
u(x) = 1.0 + 2*x[1] - x[2] + 3*x[1]^2 - x[1] * x[2] + 0.5*x[2]^2
∇u(x) = VectorValue(2 + 6*x[1] - x[2], -1 - x[1] + x[2])
@test evaluate(dofs, GenericField(u))[4:6] ≈ exact_edge_dofs(∇u)

# ... and only an approximation otherwise
w(x) = sin(2*x[1]) * cos(3*x[2])
∇w(x) = VectorValue(2*cos(2*x[1]) * cos(3*x[2]), -3*sin(2*x[1]) * sin(3*x[2]))
got = evaluate(dofs, GenericField(w))[4:6]
want = exact_edge_dofs(∇w)
@test !isapprox(got, want)
@test isapprox(got, want; rtol=0.05)

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
