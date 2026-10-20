module HermiteRefFEsTests

using LinearAlgebra
using Combinatorics: permutations
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

# Transposed Jacobian of the affine map from the reference D-simplex onto the
# simplex with vertices `verts`, in Gridap's convention Jt[i,j] = ∂Fⱼ/∂ξⁱ.
function jacobian_t(verts)
  D = length(verts) - 1
  # column-major data: column i is the edge vector verts[i+1] - verts[1]
  J = TensorValue{D,D}(ntuple(k -> (verts[(k-1)÷D+2] - verts[1])[(k-1)%D+1], D*D))
  transpose(J)
end

affine_map(verts) = ξ -> verts[1] + sum((verts[i+1] - verts[1]) * ξ[i] for i in 1:length(verts)-1)

# Generalized Vandermonde matrix Q[l,i] = ℓₗ(F*(Ψ̂ᵢ)) of [Kirby, (3.10)]: the
# physical Hermite DoFs of the simplex `verts` -- values at the vertices and the
# 2-face barycenters, then vertex gradients in the global Cartesian frame, in the
# element's own DoF order -- evaluated on the pullback of the reference shape
# functions. Computed from scratch, independently of the change-of-basis code.
function hermite_vandermonde(reffe, verts)
  p = get_polytope(reffe)
  D = num_dims(p)
  Ψ̂ = get_shapefuns(reffe)
  ndofs = num_dofs(reffe)

  K = inv(jacobian_t(verts))   # J⁻ᵀ: ∇u = K∇û
  V = eltype(evaluate(Ψ̂, [first(get_vertex_coordinates(p))]))
  comps = component_basis(V)
  nc = length(comps)

  nodes = copy(get_vertex_coordinates(p))
  D ≥ 2 && append!(nodes, [sum(c) / length(c) for c in get_face_coordinates(p, 2)])
  nn = length(nodes)
  nv = num_faces(p, 0)

  Q = zeros(Float64, ndofs, ndofs)
  vals = evaluate(Ψ̂, nodes)
  for c in 1:nc, n in 1:nn, j in 1:ndofs
    Q[(c-1)*nn+n, j] = vals[n, j] ⊙ comps[c]     # values, component by component
  end
  ∇Ψ̂ = evaluate(Broadcasting(∇)(Ψ̂), get_vertex_coordinates(p))
  off = nc * nn
  for v in 1:nv, j in 1:ndofs
    g = K ⋅ ∇Ψ̂[v, j]                              # physical gradient, g[d, c…]
    for d in 1:D, c in 1:nc
      eᵈ = VectorValue(ntuple(i -> i == d ? 1.0 : 0.0, D))
      Q[off+(v-1)*D*nc+(d-1)*nc+c, j] = (eᵈ ⋅ g) ⊙ comps[c]
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
# cycling through all the permutations of the cell vertices. The mesh is
# geometrically identical; only the local vertex numbering changes.
function permute_cells(model)
  grid = get_grid(model)
  cell_nodes = get_cell_node_ids(grid)
  perms = collect(permutations(1:num_cell_dims(model)+1))
  permuted = [collect(cell_nodes[i])[perms[mod1(i, length(perms))]]
              for i in 1:length(cell_nodes)]
  pgrid = UnstructuredGrid(get_node_coordinates(grid), Table(permuted),
                           get_reffes(grid), get_cell_type(grid))
  tag_boundary!(UnstructuredDiscreteModel(pgrid))
end

unit_segment(n) = CartesianDiscreteModel((0, 1), (n,))
unit_square(n) = simplexify(CartesianDiscreteModel((0, 1, 0, 1), (n, n)))
unit_cube(n) = simplexify(CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (n, n, n)))

T = Float64

############################################################################################
# Reference element
############################################################################################

for (p, ndofs) in ((SEGMENT, 4), (TRI, 10), (TET, 20))
  reffe = HermiteRefFE(T, p)

  @test reffe isa GenericRefFE
  @test get_name(reffe) isa Hermite
  @test num_dofs(reffe) == ndofs
  @test length(get_prebasis(reffe)) == ndofs
  @test Conformity(reffe) == H1Conformity()
  @test Pushforward(Hermite) == IdentityPiolaMap()
  @test reffe == ReferenceFE(p, hermite, T)
  @test reffe == ReferenceFE(p, hermite, T, 3)
  test_reference_fe(reffe)

  # Unisolvency
  E = evaluate(get_dof_basis(reffe), get_shapefuns(reffe))
  @test E ≈ Matrix(I, ndofs, ndofs)
  @test cond(evaluate(get_dof_basis(reffe), get_prebasis(reffe))) < 1e3
end

# The element is defined on simplices only, at order 3 only
@test_throws ErrorException HermiteRefFE(T, QUAD)
@test_throws ErrorException HermiteRefFE(T, HEX)
@test_throws ErrorException ReferenceFE(TRI, hermite, T, 2)

# One value and D gradient components per vertex, one value per 2-face
hermite1 = HermiteRefFE(T, SEGMENT)
@test get_face_own_dofs(hermite1) == [[1, 3], [2, 4], []]

hermite2 = HermiteRefFE(T, TRI)
@test get_face_own_dofs(hermite2) == [
#  φ  ∂x(φ) ∂y(φ)  at v1
  [1, 5,    6],
#  φ  ∂x(φ) ∂y(φ)  at v2
  [2, 7,    8],
#  φ  ∂x(φ) ∂y(φ)  at v3
  [3, 9,    10],
  [], [], [],
#  φ at TRI barycenter
  [4]
]
@test get_face_dofs(hermite2) == [
  [1, 5, 6], [2, 7, 8], [3, 9, 10],
  [1, 5, 6, 2, 7, 8], [1, 5, 6, 3, 9, 10], [2, 7, 8, 3, 9, 10],
  [1, 5, 6, 2, 7, 8, 3, 9, 10, 4]
]

hermite3 = HermiteRefFE(T, TET)
@test get_face_own_dofs(hermite3) == [
  [1, 9, 10, 11], [2, 12, 13, 14], [3, 15, 16, 17], [4, 18, 19, 20],
  [], [], [], [], [], [],
  [5], [6], [7], [8],           # the barycenters of the four faces
  []
]

# The shape functions in the monomial basis, against DefElement
dofs = get_dof_basis(hermite1)
P3_1D_monoms = FEEC_poly_basis(Val(1), T, 3, 0, :P, Monomial)
shapefuns_monom_coordinates = inv(evaluate(dofs, P3_1D_monoms))
@test shapefuns_monom_coordinates ≈ T[
  1  -0   0  -0;
  0   0   1  -0;
 -3   3  -2  -1;
  2  -2   1   1;
] # from https://defelement.org/elements/examples/interval-hermite-3.html

dofs = get_dof_basis(hermite2)
P3_2D_monoms = FEEC_poly_basis(Val(2), T, 3, 0, :P, Monomial)
shapefuns_monom_coordinates = inv(evaluate(dofs, P3_2D_monoms))
@test shapefuns_monom_coordinates ≈ T[
   1  -0   0    0   0   0   0  -0   0  -0;
   0   0   0    0   1   0   0  -0   0  -0;
  -3   3   0    0  -2   0  -1  -0   0  -0;
   2  -2   0    0   1   0   1  -0   0  -0;
   0   0   0    0   0   1   0  -0   0  -0;
 -13  -7  -7   27  -3  -3   2  -1  -1   2;
  13   7   7  -27   3   2  -2   2   1  -2;
  -3   0   3   -0   0  -2   0   0   0  -1;
  13   7   7  -27   2   3  -2   1   2  -2;
   2   0  -2    0   0   1   0   0   0   1;
] # from https://defelement.org/elements/examples/triangle-hermite-3.html

############################################################################################
# Cartesian product elements
############################################################################################

V = SymTracelessTensorValue{2,T}
hermite1_q = HermiteRefFE(V, SEGMENT)
dofs = get_dof_basis(hermite1_q)
shapefuns = get_shapefuns(hermite1_q)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-15
test_reference_fe(hermite1_q)

P3_1D_monoms = FEEC_poly_basis(Val(1), V, 3, 0, :P, Monomial; cart_prod=true)
shapefuns_monom_coordinates = inv(evaluate(dofs, P3_1D_monoms))
@test shapefuns_monom_coordinates ≈ T[
# v1  v2  v1  v2  v1   v1   v2   v2
# φ₁  φ₁  φ₂  φ₂ ∂xφ₁ ∂xφ₂ ∂xφ₁ ∂xφ₂
  1   0  -0  -0   0    0    0   -0;
  0   0   1  -0   0    0    0   -0;
  0   0   0  -0   1    0    0   -0;
  0   0   0   0   0    1    0   -0;
 -3   3   0   0  -2   -0   -1   -0;
  0   0  -3   3   0   -2    0   -1;
  2  -2   0   0   1    0    1   -0;
  0   0   2  -2   0    1    0    1;
]

@test get_face_own_dofs(hermite1_q) == Vector{Int}[
#  φ₁ φ₂ ∂x(φ₁) ∂x(φ₂)  at v1
  [1, 3,  5,    6],
#  φ₁ φ₂ ∂x(φ₁) ∂x(φ₂)  at v2
  [2, 4,  7,    8],
  []
]

V = VectorValue{2,T}
hermite2_v = HermiteRefFE(V, TRI)
@test num_dofs(hermite2_v) == 20
@test norm(evaluate(get_dof_basis(hermite2_v), get_shapefuns(hermite2_v)) - I) <= 1.e-15
test_reference_fe(hermite2_v)
@test get_face_own_dofs(hermite2_v) == Vector{Int}[
#  φ₁ φ₂ ∂x(φ₁) ∂x(φ₂) ∂y(φ₁) ∂y(φ₂)  at v1
  [1, 5,  9,    10,    11,    12],
  [2, 6, 13,    14,    15,    16],
  [3, 7, 17,    18,    19,    20],
  [], [], [],
#  φ₁ φ₂  at the barycenter
  [4, 8]
]

hermite3_v = HermiteRefFE(V, TET)
dofs = get_dof_basis(hermite3_v)
shapefuns = get_shapefuns(hermite3_v)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-15
test_reference_fe(hermite3_v)
@test get_face_own_dofs(hermite3_v) == Vector{Int}[
#  φ₁ φ₂ ∂x(φ₁) ∂x(φ₂) ∂y(φ₁) ∂y(φ₂) ∂z(φ₁) ∂z(φ₂)  at v1
  [1, 9,  17, 18, 19, 20, 21, 22],
  [2, 10, 23, 24, 25, 26, 27, 28],
  [3, 11, 29, 30, 31, 32, 33, 34],
  [4, 12, 35, 36, 37, 38, 39, 40],
  [], [], [], [], [], [],
#  φ₁ φ₂  at face 1 barycenter
  [5, 13],
#  φ₁ φ₂  at face 2 barycenter
  [6, 14],
#  φ₁ φ₂  at face 3 barycenter
  [7, 15],
#  φ₁ φ₂  at face 4 barycenter
  [8, 16],
  []
]

############################################################################################
# The DoFs are the claimed functionals
############################################################################################

# Scalar: vertex values, the barycenter value, then the gradient components at
# each vertex
prebasis = get_prebasis(hermite2)
Vd = evaluate(get_dof_basis(hermite2), prebasis)
V̂ = get_vertex_coordinates(TRI)
b̂ = sum(V̂) / 3
@test Vd[1:3, :] ≈ evaluate(prebasis, V̂)
@test Vd[4:4, :] ≈ evaluate(prebasis, [b̂])
∇φ = evaluate(Broadcasting(∇)(prebasis), V̂)
for v in 1:3
  @test Vd[4+2*v-1, :] ≈ [g[1] for g in view(∇φ, v, :)]
  @test Vd[4+2*v, :] ≈ [g[2] for g in view(∇φ, v, :)]
end

# Vector-valued: the same, component after component for the values, and
# direction by direction with the component running fastest for the gradients
prebasis = get_prebasis(hermite2_v)
Vd = evaluate(get_dof_basis(hermite2_v), prebasis)
φv = evaluate(prebasis, [V̂..., b̂])
for c in 1:2
  @test Vd[4*(c-1) .+ (1:4), :] ≈ [φ[c] for φ in φv]
end
∇φ = evaluate(Broadcasting(∇)(prebasis), V̂)
for v in 1:3, d in 1:2, c in 1:2
  @test Vd[8+4*(v-1)+2*(d-1)+c, :] ≈ [g[d, c] for g in view(∇φ, v, :)]
end

# A linear field, whose gradient DoFs are its coefficients
u = GenericField(x -> VectorValue(x[1] + 2*x[2], 3*x[1] + 4*x[2]))
@test evaluate(get_dof_basis(hermite2_v), u)[9:12] ≈ [1.0, 3.0, 2.0, 4.0]

############################################################################################
# Transformation to a physical cell
############################################################################################

test_cells = Dict(
  SEGMENT => (
    [Point(0.0), Point(1.0)],                             # the reference cell
    [Point(0.3), Point(1.7)],
    [Point(0.3), Point(-1.2)],                            # det J < 0
  ),
  TRI => (
    [Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0)],  # the reference cell
    [Point(0.3, 0.1), Point(1.7, 0.2), Point(0.6, 2.3)],  # det J > 0
    [Point(0.0, 0.0), Point(0.0, 1.0), Point(1.0, 0.0)],  # det J < 0
    [Point(-1.0, 0.5), Point(2.0, -0.4), Point(0.1, 0.3)],
  ),
  TET => (
    [Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0)],
    [Point(0.3, 0.1, -0.2), Point(1.7, 0.2, 0.4), Point(0.6, 2.3, 0.1), Point(-0.5, 0.7, 1.9)],
    [Point(0.0, 0.0, 0.0), Point(0.0, 1.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 0.0, 1.0)],
  ),
)

# The vector-valued cells have a non-symmetric Jt, which is what tells the
# gradient block Kg ⊗ I from I ⊗ Kg and Kg from Kgᵀ.
for p in (SEGMENT, TRI, TET), V in (T, VectorValue{2,T})
  reffe = HermiteRefFE(V, p)
  n = num_dofs(reffe)
  for (i, verts) in enumerate(test_cells[p])
    Jt = jacobian_t(verts)
    P = copy(evaluate(ReferenceFEs.HermiteChangeOfBasis(reffe, false), Jt))
    Pinvt = copy(evaluate(ReferenceFEs.HermiteChangeOfBasis(reffe, true), Jt))

    @test transpose(Pinvt) * P ≈ Matrix(I, n, n)

    Q = hermite_vandermonde(reffe, verts)
    @test Q * P ≈ Matrix(I, n, n)

    # the values are preserved, and the gradients are not unless J = I
    nvalues = n - num_faces(p, 0) * num_dims(p) * num_indep_components(V)
    @test P[1:nvalues, :] ≈ Matrix(I, n, n)[1:nvalues, :]
    if i == 1
      @test P ≈ Matrix(I, n, n)
    else
      @test !isapprox(P, Matrix(I, n, n))
    end
  end
end

# Against the element built directly on the physical cell: its shape functions
# at F(x̂) are the transformed reference ones at x̂.
for p in (SEGMENT, TRI, TET), V in (T, VectorValue{2,T})
  reffe = HermiteRefFE(V, p)
  x̂ = [Point(ntuple(i -> 0.1 * i + 0.05, num_dims(p))), Point(ntuple(i -> 0.2 / i, num_dims(p)))]
  for verts in test_cells[p]
    phys = HermiteRefFE(V, p; vertices=verts)
    @test get_polytope(phys) == p
    @test get_face_own_dofs(phys) == get_face_own_dofs(reffe)
    @test norm(evaluate(get_dof_basis(phys), get_shapefuns(phys)) - I) < 1e-12

    P = copy(evaluate(ReferenceFEs.HermiteChangeOfBasis(reffe, false), jacobian_t(verts)))
    Ψ = evaluate(linear_combination(P, get_shapefuns(reffe)), x̂)
    Ψphys = evaluate(get_shapefuns(phys), affine_map(verts).(x̂))
    @test maximum(norm, Ψ - Ψphys) < 1e-12   # many entries are exact zeros
  end
end

############################################################################################
# FE space
#
# The discrete space must not depend on how each cell happens to list its
# vertices, so this runs on a mesh whose cells are sorted by global vertex id and
# on the same mesh with the cell vertices permuted.
############################################################################################

function test_hermite_fe_space(model, ::Type{V}) where V
  D = num_cell_dims(model)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 8)
  Vh = FESpace(model, ReferenceFE(hermite, V, 3))

  topo = get_grid_topology(model)
  nc = num_indep_components(V)
  n2faces = D ≥ 2 ? num_faces(topo, 2) : 0
  @test num_free_dofs(Vh) == nc * ((D + 1) * num_faces(topo, 0) + n2faces)

  # The space is P₃ on each cell and the DoFs are unisolvent, so a global cubic
  # must be interpolated exactly.
  s(x) = 1 + x[1] - 2*x[min(2, D)] + x[1]^2 * x[min(2, D)] - 0.5*x[1]^3 + x[min(3, D)]^3
  u = V <: MultiValue ? (x -> VectorValue(s(x), 2*s(x))) : s
  uh = interpolate(u, Vh)
  eu = uh - u
  @test sqrt(sum(∫(eu ⊙ eu)dΩ)) < 1e-13

  # A field outside the space, to test conformity
  t(x) = sin(2*x[1]) * cos(3*x[min(2, D)])
  w = V <: MultiValue ? (x -> VectorValue(t(x), 2*t(x))) : t
  wh = interpolate(w, Vh)
  ew = wh - w
  @test sqrt(sum(∫(ew ⊙ ew)dΩ)) > 1e-4

  # Hermite is C⁰ ...
  if D ≥ 2
    Λ = SkeletonTriangulation(model)
    dΛ = Measure(Λ, 8)
    @test sqrt(sum(∫(jump(wh) ⊙ jump(wh))dΛ)) < 1e-12
    # ... but not C¹ along the facets ...
    @test sqrt(sum(∫(jump(∇(wh)) ⊙ jump(∇(wh)))dΛ)) > 1e-2
  end

  # ... only at the vertices: the gradient of `wh` at a vertex is the same from
  # every cell touching it.
  p = get_polytope(only(get_reffes(model)))
  cell_vertices = Geometry.get_faces(topo, D, 0)
  ξ = get_vertex_coordinates(p)
  cell_∇wh = lazy_map(evaluate, get_data(∇(wh)), fill(ξ, num_cells(model)))
  vertex_grads = Dict{Int,Any}()
  for cell in 1:num_cells(model)
    g = cell_∇wh[cell]
    for (lv, v) in enumerate(cell_vertices[cell])
      gv = get!(vertex_grads, v, g[lv])
      @test norm(g[lv] - gv) < 1e-12
    end
  end
end

test_hermite_fe_space(unit_segment(4), T)
test_hermite_fe_space(unit_segment(4), VectorValue{2,T})
test_hermite_fe_space(unit_square(3), T)
test_hermite_fe_space(permute_cells(unit_square(3)), T)
test_hermite_fe_space(unit_square(3), VectorValue{2,T})
test_hermite_fe_space(permute_cells(unit_square(3)), VectorValue{2,T})
test_hermite_fe_space(unit_cube(2), T)
test_hermite_fe_space(permute_cells(unit_cube(2)), T)

# The permuted mesh does permute cell vertices
@test all(issorted, Geometry.get_faces(get_grid_topology(unit_square(3)), 2, 0))
@test !all(issorted, Geometry.get_faces(get_grid_topology(permute_cells(unit_square(3))), 2, 0))
@test !all(issorted, Geometry.get_faces(get_grid_topology(permute_cells(unit_cube(2))), 3, 0))

# A Poisson solve, whose solution must not depend on the cell vertex order either
function solve_poisson(model)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 8)
  uex(x) = sin(2π*x[1]) * sin(2π*x[2])
  f(x) = 8π^2 * uex(x)
  Vh = FESpace(model, ReferenceFE(hermite, T, 3); dirichlet_tags="boundary")
  U = TrialFESpace(Vh, uex)
  a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
  l(v) = ∫(f * v)dΩ
  uh = solve(AffineFEOperator(a, l, U, Vh))
  e = uh - uex
  sqrt(sum(∫(e * e)dΩ))
end
e_sorted = solve_poisson(unit_square(8))
e_permuted = solve_poisson(permute_cells(unit_square(8)))
@test e_sorted < 2e-3
@test e_sorted ≈ e_permuted

end # module
