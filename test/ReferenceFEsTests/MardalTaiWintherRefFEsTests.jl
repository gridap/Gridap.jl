module MardalTaiWintherRefFEsTests

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
# Helpers -- meshes, shared
############################################################################################

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

# Rebuild `model` with the vertices of each cell listed in a permuted order. The
# mesh is geometrically identical; only the local vertex numbering -- and hence
# the direction in which each cell traverses its edges (2D) or the order in which
# it lists each face's vertices (3D) -- changes.
function permute_cells(model, perms)
  grid = get_grid(model)
  cell_nodes = get_cell_node_ids(grid)
  permuted = [collect(cell_nodes[i])[perms[mod1(i, length(perms))]]
              for i in 1:length(cell_nodes)]
  pgrid = UnstructuredGrid(get_node_coordinates(grid), Table(permuted),
                           get_reffes(grid), get_cell_type(grid))
  tag_boundary!(UnstructuredDiscreteModel(pgrid))
end

# even and odd permutations, so both orientations occur
permute_cells_2d(model) = permute_cells(model, ([1,2,3],[2,3,1],[3,1,2],
                                                [2,1,3],[1,3,2],[3,2,1]))
permute_cells_3d(model) = permute_cells(model, ([1,2,3,4],[2,3,4,1],[4,1,2,3],
                                                [2,1,3,4],[1,3,2,4],[4,3,2,1]))

unit_square(n) = simplexify(CartesianDiscreteModel((0, 1, 0, 1), (n, n)))
unit_cube(n) = simplexify(CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (n, n, n)))

############################################################################################
# Helpers -- 2D
############################################################################################

# Quadrature data for integrating over edge `k` of `p`: the points mapped to the
# cell, the weights scaled by the edge length (so they integrate `ds`), the
# coordinates `ŝ ∈ [0,1]` on the reference segment, and the edge frame.
function edge_quadrature(p::Polytope{2}, k::Integer, degree::Integer)
  quad = Quadrature(SEGMENT, degree)
  ŝ = get_coordinates(quad)
  ŵ = get_weights(quad)

  v1, v2 = get_face_coordinates(p, 1)[k]
  e = v2 - v1
  L = norm(e)
  t = e / L

  x = [v1 + si[1] * e for si in ŝ]
  (x, ŵ .* L, ŝ, t, ReferenceFEs._rot90(t))
end

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

# Generalized Vandermonde matrix Q[l,i] = ℓₗ(F*(Ψ̂ᵢ)) of
# [Aznaran, Farrell & Kirby, (5.10)]: the physical 2D MTW DoFs of the triangle
# `verts`, evaluated on the contravariant Piola push-forward of the reference
# shape functions. Computed from scratch by quadrature on the physical edges,
# independently of the change-of-basis code.
function mtw_vandermonde(reffe, verts; degree=8)
  p = get_polytope(reffe)
  Ψ̂ = get_shapefuns(reffe)
  ndofs = num_dofs(reffe)

  Jt = jacobian_t(verts)
  detJ = det(Jt)
  F = affine_map(verts)

  quad = Quadrature(SEGMENT, degree)
  ŝ = get_coordinates(quad)
  ŵ = get_weights(quad)
  μ = evaluate(LegendreBasis(Val(1), Float64, 1), ŝ)

  Q = zeros(Float64, ndofs, ndofs)
  for k in 1:num_faces(p, 1)
    v̂a, v̂b = get_face_coordinates(p, 1)[k]
    xa, xb = F(v̂a), F(v̂b)
    L = norm(xb - xa)
    t = (xb - xa) / L
    n = ReferenceFEs._rot90(t)

    x̂ = [v̂a + si[1] * (v̂b - v̂a) for si in ŝ]
    Ψ̂x = evaluate(Ψ̂, x̂)                                  # (nq, ndofs)
    Φx = map(v -> (1 / abs(detJ)) * (v ⋅ Jt), Ψ̂x)         # push-forward

    o = 3*(k - 1)
    for j in 1:ndofs
      Φn = map(φ -> φ ⋅ n, view(Φx, :, j))
      Φt = map(φ -> φ ⋅ t, view(Φx, :, j))
      Q[o+1, j] = sum(ŵ .* L .* Φn .* μ[:, 1])
      Q[o+2, j] = sum(ŵ .* L .* Φn .* μ[:, 2])
      Q[o+3, j] = sum(ŵ .* L .* Φt .* μ[:, 1])
    end
  end
  Q
end

############################################################################################
# Helpers -- 3D
############################################################################################

vec3(v) = [v[1], v[2], v[3]]

# Geometry of a face as a triangle in R³.
function face_frame(â, b̂, ĉ)
  ê1, ê2 = b̂ - â, ĉ - â
  nvec = cross(ê1, ê2)
  area2 = norm(nvec)
  n = nvec / area2
  t1 = ê1 / norm(ê1)
  (ê1, ê2, area2, n, t1, cross(n, t1), (â + b̂ + ĉ) / 3)
end

# The 24 generators of V(T) = P₁(T;R³) + curl(b P₁(T;R³)), written the direct
# way: the quartic bubble as a closure, the twelve vector linears, and the curl
# of the bubble times a vector linear taken by automatic differentiation. The
# element writes the Bernstein coefficients down exactly instead, so these are
# the test's independent oracle -- they share no code with the construction they
# check, and AD through a closure is unobjectionable in a test.
function tw_bubble(p::Polytope{3})
  @assert get_vertex_coordinates(p) == [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                                        Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0)]
  x -> (1 - x[1] - x[2] - x[3]) * x[1] * x[2] * x[3]
end

function tw_generators(p::Polytope{3})
  b = tw_bubble(p)
  lins = (x -> 1.0, x -> x[1], x -> x[2], x -> x[3])
  es = (VectorValue(1.0, 0.0, 0.0), VectorValue(0.0, 1.0, 0.0), VectorValue(0.0, 0.0, 1.0))
  linear = [(x -> l(x) * e) for l in lins, e in es][:]
  bubbly = [(x -> b(x) * l(x) * e) for l in lins, e in es][:]
  vcat(linear, [(x -> (cross(∇, f))(x)) for f in bubbly])
end

# The full 24x24 generalized Vandermonde of the physical 3D DoFs of the
# tetrahedron `X`, computed from scratch: each face is integrated with its
# vertices in the order given by `pids`, matching what the change of basis uses.
# Computed over *all* 24 columns, so it also tests that Q is block diagonal.
function tw_vandermonde(reffe, X, pids; degree=8)
  p = get_polytope(reffe)
  Ψ = get_shapefuns(reffe)
  own = get_face_own_dofs(reffe)
  fdim = get_dimrange(p, 2)
  fdofs = [own[f] for f in fdim]
  fverts = get_faces(p, 2, 0)
  vc = get_vertex_coordinates(p)

  Jm = hcat(vec3(X[2] - X[1]), vec3(X[3] - X[1]), vec3(X[4] - X[1]))
  dJ = det(Jm)
  quad = Quadrature(TRI, degree)
  qc, qw = get_coordinates(quad), get_weights(quad)

  Q = zeros(24, 24)
  for (fi, dofs) in enumerate(fdofs)
    lv = fverts[fi][invperm(get_face_vertex_permutations(p, 2)[fi][pids[fi]])]
    â, b̂, ĉ = vc[lv[1]], vc[lv[2]], vc[lv[3]]
    x̂ = [â + q[1] * (b̂ - â) + q[2] * (ĉ - â) for q in qc]
    Ψx = evaluate(Ψ, x̂)

    # the element's weights: scaled normal n = e₁×e₂ and the tangential triple
    # n×e₁, n×e₂, n×(x-a), integrated over the parameter domain
    av, bv, cv = X[lv[1]], X[lv[2]], X[lv[3]]
    e1, e2 = bv - av, cv - av
    n = cross(e1, e2)
    w4, w5 = cross(n, e1), cross(n, e2)
    for (q, wq) in enumerate(qw)
      u, v = qc[q][1], qc[q][2]
      w6 = cross(n, u * e1 + v * e2)
      for i in 1:24
        φ = VectorValue(Jm * vec3(Ψx[q, i]) / abs(dJ))
        φn = φ ⋅ n
        Q[dofs[1], i] += wq * φn
        Q[dofs[2], i] += wq * φn * u
        Q[dofs[3], i] += wq * φn * v
        Q[dofs[4], i] += wq * (φ ⋅ w4)
        Q[dofs[5], i] += wq * (φ ⋅ w5)
        Q[dofs[6], i] += wq * (φ ⋅ w6)
      end
    end
  end
  Q, fdofs
end

############################################################################################
# 2D -- reference element
############################################################################################

reffe2 = MardalTaiWintherRefFE(Float64, TRI)

@test reffe2 isa GenericRefFE
@test get_name(reffe2) isa MardalTaiWinther
@test num_dofs(reffe2) == 9
# the prebasis is the ambient P₃(K;R²), not the element's own 9-dimensional
# space: the augmented construction restricts to the DoFs through the shape
# functions, so `num_dofs != length(prebasis)` here by design
@test length(get_prebasis(reffe2)) == 20
@test Conformity(reffe2) == DivConformity()
@test Pushforward(MardalTaiWinther) == ContraVariantPiolaMap()
@test reffe2 == ReferenceFE(TRI, mtw, Float64)
test_reference_fe(reffe2)

# Three DoFs per edge, none on the vertices or the cell
@test length.(get_face_own_dofs(reffe2)) == [0, 0, 0, 3, 3, 3, 0]

# Duality of the DoF and shape function bases
@test evaluate(get_dof_basis(reffe2), get_shapefuns(reffe2)) ≈ Matrix(I, 9, 9)

@test_throws ErrorException MardalTaiWintherRefFE(Float64, QUAD)
@test_throws ErrorException MardalTaiWintherRefFE(Float64, HEX)
@test_throws ErrorException ReferenceFE(TRI, mtw, Float64, 2)

############################################################################################
# 2D -- the space is what the constraints say it is
############################################################################################

Ψ2 = get_shapefuns(reffe2)

# div Ψ ∈ P₀(K)
divΨ = map(tr, evaluate(Broadcasting(∇)(Ψ2), get_coordinates(Quadrature(TRI, 6))))
for j in 1:9
  @test all(d -> isapprox(d, divΨ[1, j]; atol=1e-11), view(divΨ, :, j))
end

# (Ψ⋅n)|ₑ ∈ P₁(e): the degree 2 and 3 Legendre moments vanish
μbasis = LegendreBasis(Val(1), Float64, 3)
for k in 1:3
  xe, we, ŝ, _, n = edge_quadrature(TRI, k, 8)
  μ = evaluate(μbasis, ŝ)
  Ψn = map(φ -> φ ⋅ n, evaluate(Ψ2, xe))
  for j in 1:9, i in 3:4
    @test abs(sum(we .* Ψn[:, j] .* μ[:, i])) < 1e-12
  end
end

# P₁(K;R²) ⊂ MTW(K)
xs2 = get_coordinates(Quadrature(TRI, 4))
Ψx2 = evaluate(Ψ2, xs2)
A2 = vcat([Ψx2[i, j][1] for i in axes(Ψx2, 1), j in axes(Ψx2, 2)],
          [Ψx2[i, j][2] for i in axes(Ψx2, 1), j in axes(Ψx2, 2)])
for u in (x -> VectorValue(1.0, 0.0), x -> VectorValue(0.0, 1.0),
          x -> VectorValue(x[1], 0.0), x -> VectorValue(x[2], 0.0),
          x -> VectorValue(0.0, x[1]), x -> VectorValue(0.0, x[2]))
  b = vcat([u(x)[1] for x in xs2], [u(x)[2] for x in xs2])
  @test norm(A2 * (A2 \ b) - b) < 1e-11
end

############################################################################################
# 2D -- transformation to a physical cell
############################################################################################

ts2, ns2 = ReferenceFEs._edge_frames(TRI)

test_cells_2d = (
  [Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0)],  # the reference cell
  [Point(0.3, 0.1), Point(1.7, 0.2), Point(0.6, 2.3)],  # det J > 0
  [Point(0.0, 0.0), Point(0.0, 1.0), Point(1.0, 0.0)],  # det J < 0
  [Point(-1.0, 0.5), Point(2.0, -0.4), Point(0.1, 0.3)],
)
# every way the cell's edge traversals can disagree with the global ones
test_signs_2d = ((1.0, 1.0, 1.0), (-1.0, 1.0, 1.0), (1.0, -1.0, 1.0),
                 (-1.0, -1.0, -1.0))

for verts in test_cells_2d, σ in test_signs_2d
  Jt = jacobian_t(verts)
  P = copy(evaluate(FESpaces.MTWChangeOfBasis(ts2, ns2, false), Jt, σ))
  Pinvt = copy(evaluate(FESpaces.MTWChangeOfBasis(ts2, ns2, true), Jt, σ))

  # The closed forms are consistent: Pinvt is the transposed inverse of P
  @test transpose(Pinvt) * P ≈ Matrix(I, 9, 9)

  # The physical nodal basis is Ψ = linear_combination(P, F*(Ψ̂)), i.e. the
  # physical DoFs are dual to it: Q P = I with Q the Vandermonde matrix.
  # Reversing the global direction of an edge flips its ℓⁿ⁰ and ℓᵗ⁰ rows.
  Dσ = Diagonal(vcat(([σe, 1.0, σe] for σe in σ)...))
  Q = Dσ * mtw_vandermonde(reffe2, verts)
  @test Q * P ≈ Matrix(I, 9, 9)

  # Without the change of basis, the pushed-forward basis is NOT nodal
  if verts != test_cells_2d[1] || σ != test_signs_2d[1]
    @test !isapprox(Q, Matrix(I, 9, 9))
  end
end

############################################################################################
# 2D -- FE space
#
# The discrete space must not depend on how each cell happens to list its
# vertices, so this runs on a sorted mesh and on the same mesh permuted.
############################################################################################

function test_mtw_fe_space_2d(model)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 8)
  V = FESpace(model, ReferenceFE(mtw, Float64, 1))

  topo = get_grid_topology(model)
  @test num_free_dofs(V) == 3*num_faces(topo, 1)

  # MTW contains P₁(K;R²), so a global linear field is interpolated exactly
  u(x) = VectorValue(1.0 + 2.0*x[1] - x[2], 3.0 - x[1] + 4.0*x[2])
  uh = interpolate(u, V)
  @test sqrt(sum(∫((uh - u) ⋅ (uh - u))dΩ)) < 1e-11

  # A field outside the space, to test conformity
  w(x) = VectorValue(sin(3*x[1]) * x[2]^2, cos(2*x[2]) + x[1]^3)
  wh = interpolate(w, V)
  @test sqrt(sum(∫((wh - w) ⋅ (wh - w))dΩ)) > 1e-6

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 8)
  n_Λ = get_normal_vector(Λ)
  t_Λ = Operation(ReferenceFEs._rot90)(n_Λ)

  # H(div)-conforming: the normal component is continuous
  @test sqrt(sum(∫(jump(wh ⋅ n_Λ) * jump(wh ⋅ n_Λ))dΛ)) < 1e-10

  # H¹-nonconforming: only the mean tangential trace is continuous
  @test maximum(abs, get_array(∫(jump(wh ⋅ t_Λ))dΛ)) < 1e-11
  @test sqrt(sum(∫(jump(wh ⋅ t_Λ) * jump(wh ⋅ t_Λ))dΛ)) > 1e-6

  # and the interpolant is the same function on both meshes
  @test sqrt(sum(∫((wh - w) ⋅ (wh - w))dΩ)) ≈ 0.04430864610837 atol=1e-12
end

test_mtw_fe_space_2d(unit_square(3))
test_mtw_fe_space_2d(permute_cells_2d(unit_square(3)))

############################################################################################
# 2D -- edge orientation
############################################################################################

# `simplexify(CartesianDiscreteModel(...))` already lists every cell's vertices
# in increasing global id order, so no edge needs flipping
sorted_model = unit_square(3)
@test all(issorted, Geometry.get_faces(get_grid_topology(sorted_model), 2, 0))
@test all(σ -> all(isequal(1.0), σ), FESpaces._edge_signs(sorted_model, TRI))

# After permuting, a sign is -1 exactly where the cell lists the edge against the
# order the mesh stores it in. The master is the mesh's own order, not `sortperm`
# of the global ids: the two agree on a sorted mesh but not after permuting, and
# differ only by a global sign per edge.
pmodel = permute_cells_2d(sorted_model)
ptopo = get_grid_topology(pmodel)
cell_vertices = Geometry.get_faces(ptopo, 2, 0)
@test !all(issorted, cell_vertices)

signs = FESpaces._edge_signs(pmodel, TRI)
edge_lvertices = get_faces(TRI, 1, 0)
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

############################################################################################
# 3D -- reference element
############################################################################################

reffe3 = MardalTaiWintherRefFE(Float64, TET)

@test reffe3 isa GenericRefFE
@test get_name(reffe3) isa MardalTaiWinther
@test num_dofs(reffe3) == 24
# in 3D the space is a sum, not a kernel, so the prebasis *is* the element's own
# space and no constraints are needed
@test length(get_prebasis(reffe3)) == 24
@test Conformity(reffe3) == DivConformity()
@test reffe3 == ReferenceFE(TET, mtw, Float64)
test_reference_fe(reffe3)

# six DoFs on each of the four faces, none anywhere else
@test length.(get_face_own_dofs(reffe3)) == [zeros(Int, 10); fill(6, 4); 0]

@test evaluate(get_dof_basis(reffe3), get_shapefuns(reffe3)) ≈ Matrix(I, 24, 24)

@test_throws ErrorException ReferenceFE(TET, mtw, Float64, 2)

############################################################################################
# 3D -- the space is the claimed sum
############################################################################################

Ψ3 = get_shapefuns(reffe3)
xs3 = get_coordinates(Quadrature(TET, 8))

# the prebasis really is P₁(T;R³) + curl(b P₁(T;R³)), of dimension 24
G3 = hcat([vcat([vec3(g(x)) for x in xs3]...) for g in tw_generators(TET)]...)
@test rank(G3) == 24

# ... and the shape functions span the same space
Ψx3 = evaluate(Ψ3, xs3)
S3 = hcat([vcat([vec3(Ψx3[i, j]) for i in eachindex(xs3)]...) for j in 1:24]...)
@test rank(hcat(G3, S3)) == 24

# P₁(T;R³) ⊂ V(T)
for u in (x -> VectorValue(1.0, 0.0, 0.0), x -> VectorValue(0.0, x[3], 0.0),
          x -> VectorValue(x[1], 0.0, x[2]))
  b = vcat([vec3(u(x)) for x in xs3]...)
  @test norm(S3 * (S3 \ b) - b) < 1e-10
end

# the normal trace is linear on each face: its residual against span{1,u,v}
# vanishes
quad_f = Quadrature(TRI, 8)
qc, qw = get_coordinates(quad_f), get_weights(quad_f)
vc3 = get_vertex_coordinates(TET)
for (fi, lv) in enumerate(get_faces(TET, 2, 0))
  â, b̂, ĉ = vc3[lv[1]], vc3[lv[2]], vc3[lv[3]]
  ê1, ê2, area2, n, _, _, _ = face_frame(â, b̂, ĉ)
  x̂ = [â + q[1] * ê1 + q[2] * ê2 for q in qc]
  Ψf = evaluate(Ψ3, x̂)
  lin = hcat(ones(length(qc)), [q[1] for q in qc], [q[2] for q in qc])
  for j in 1:24
    vals = [Ψf[q, j] ⋅ n for q in eachindex(qc)]
    @test maximum(abs, vals - lin * (lin \ vals)) < 1e-12
  end
end

############################################################################################
# 3D -- transformation to a physical cell
############################################################################################

test_cells_3d = (
  [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0)],
  [Point(0.1, 0.2, 0.3), Point(1.3, 0.1, 0.2), Point(0.2, 1.7, 0.4), Point(0.3, 0.2, 2.1)],
  [Point(0.0, 0.0, 0.0), Point(0.0, 1.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 0.0, 1.0)],
  [Point(-1.0, 0.5, 0.2), Point(2.0, -0.4, 0.1), Point(0.1, 0.3, 1.4), Point(0.5, 2.2, -0.3)],
)
test_pids_3d = ((1, 1, 1, 1), (2, 1, 4, 1), (6, 5, 3, 2))

for X in test_cells_3d, pids in test_pids_3d
  e1, e2, e3 = X[2] - X[1], X[3] - X[1], X[4] - X[1]
  J = TensorValue(e1[1], e1[2], e1[3], e2[1], e2[2], e2[3], e3[1], e3[2], e3[3])
  Jt = transpose(J)     # columns of J are the edge vectors, so Jt is ∂Fⱼ/∂ξⁱ

  P = copy(evaluate(FESpaces.TWChangeOfBasis(reffe3, false), Jt, pids))
  Pinvt = copy(evaluate(FESpaces.TWChangeOfBasis(reffe3, true), Jt, pids))
  @test transpose(Pinvt) * P ≈ Matrix(I, 24, 24)

  Q, fdofs = tw_vandermonde(reffe3, X, pids)

  # Q is block diagonal by face -- the property the whole construction rests on
  off = maximum(maximum(abs, Q[fdofs[i], fdofs[j]]) for i in 1:4, j in 1:4 if i != j)
  @test off < 1e-10

  # and P is its inverse
  @test isapprox(Q * P, Matrix(I, 24, 24); atol=1e-10)

  if X != test_cells_3d[1]
    @test !isapprox(Q, Matrix(I, 24, 24))
  end
end

############################################################################################
# 3D -- FE space
############################################################################################

function test_mtw_fe_space_3d(model)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 8)
  V = FESpace(model, ReferenceFE(TET, mtw, Float64))

  topo = get_grid_topology(model)
  @test num_free_dofs(V) == 6*num_faces(topo, 2)

  # P₁(T;R³) ⊂ V(T), so a global linear field is interpolated exactly
  u(x) = VectorValue(1.0 + 2*x[1] - x[2], 3.0 - x[1] + x[3], 2.0 + x[2] - 3*x[3])
  uh = interpolate(u, V)
  @test sqrt(sum(∫((uh - u) ⋅ (uh - u))dΩ)) < 1e-12

  w(x) = VectorValue(sin(2*x[1]) * x[2], cos(3*x[3]) + x[1]^2, sin(x[1] + x[2] - x[3]))
  wh = interpolate(w, V)
  @test sqrt(sum(∫((wh - w) ⋅ (wh - w))dΩ)) > 1e-6

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 8)
  n_Λ = get_normal_vector(Λ)

  # H(div)-conforming: the normal component is continuous
  @test sqrt(sum(∫(jump(wh ⋅ n_Λ) * jump(wh ⋅ n_Λ))dΛ)) < 1e-12

  # nonconforming in H¹: the jump is nonzero ...
  @test sqrt(sum(∫(jump(wh) ⋅ jump(wh))dΛ)) > 1e-3

  # ... but all three rigid-motion moments of it vanish on every facet, which is
  # the weak continuity. The normal component of a constant vector pairs with
  # [u]⋅n = 0, so testing the three Cartesian directions covers the two tangential
  # translations; and since ∫_F [u] = 0, the rotation moment ∫_F [u]⋅((x-x₀)×n)
  # reduces to n⋅∫_F ([u]×x).
  for c in (VectorValue(1.0,0.0,0.0), VectorValue(0.0,1.0,0.0), VectorValue(0.0,0.0,1.0))
    @test maximum(abs, get_array(∫( jump(wh) ⋅ c )dΛ)) < 1e-12
  end
  xΛ = CellField(x -> x, Λ)
  @test maximum(abs, get_array(∫( n_Λ.⁺ ⋅ cross(jump(wh), xΛ) )dΛ)) < 1e-12
end

test_mtw_fe_space_3d(unit_cube(2))
test_mtw_fe_space_3d(permute_cells_3d(unit_cube(2)))

############################################################################################
# 3D -- face permutation ids
#
# `pindex` is Gridap's own: the permutation taking a cell's local view of a face
# to the mesh's stored vertex order for it, so `invperm` walks the face in the
# mesh's order. This checks exactly that defining property -- it is what pins the
# direction, and getting it backwards would be wrong only on the two 3-cycles.
############################################################################################

function test_mtw_face_pids(model)
  topo = get_grid_topology(model)
  cell_vertices = Geometry.get_faces(topo, 3, 0)
  face_vertices = Geometry.get_faces(topo, 2, 0)
  cell_faces = Geometry.get_faces(topo, 3, 2)
  fverts = get_faces(TET, 2, 0)
  vperms = get_face_vertex_permutations(TET, 2)

  pids = get_cell_permutations(topo, 2)
  for cell in 1:length(cell_vertices)
    v = cell_vertices[cell]
    for f in 1:4
      perm = vperms[f][pids[cell][f]]
      @test v[fverts[f][invperm(perm)]] == collect(face_vertices[cell_faces[cell][f]])
    end
  end
end

test_mtw_face_pids(unit_cube(2))
test_mtw_face_pids(permute_cells_3d(unit_cube(2)))

# a permuted mesh really does exercise the non-identity permutations
pids3 = get_cell_permutations(get_grid_topology(permute_cells_3d(unit_cube(2))), 2)
@test any(pid -> any(!isequal(1), pid), pids3)
@test all(pid -> all(p -> 1 <= p <= 6, pid), pids3)

end # module
