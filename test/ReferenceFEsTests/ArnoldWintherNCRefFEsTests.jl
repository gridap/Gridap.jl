module ArnoldWintherNCRefFEsTests

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

# Generalized Vandermonde matrix Q[l,i] = ℓₗ(F*(Ψ̂ᵢ)): the physical AWnc DoFs of
# the triangle `verts` -- the normal-normal and normal-tangential edge moments --
# evaluated on the double contravariant Piola push-forward of the reference shape
# functions, computed from scratch. The interior DoFs are *defined* as the
# push-forward of the reference ones, so their rows are the identity.
function awnc_vandermonde(reffe, verts, σ; degree=12)
  p = get_polytope(reffe)
  Ψ̂ = get_shapefuns(reffe)
  ndofs = num_dofs(reffe)
  Jt = jacobian_t(verts)
  J = transpose(Jt)
  detJ = det(Jt)
  F = affine_map(verts)
  own = get_face_own_dofs(reffe)
  nv = num_faces(p, 0)
  μb = LegendreBasis(Val(1), Float64, 1)

  Q = Matrix{Float64}(I, ndofs, ndofs)
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
        Φ = [(1 / detJ^2) * (J ⋅ ψ ⋅ transpose(J)) for ψ in view(Ψ̂x, :, j)]
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

reffe = ArnoldWintherNCRefFE(Float64, TRI)

@test reffe isa GenericRefFE
@test get_name(reffe) isa ArnoldWintherNC
@test num_dofs(reffe) == 15
# the prebasis is the ambient P₂(K;S), not the element's own 15-dimensional
# space: the augmented construction restricts to the DoFs through the shape
# functions, so `num_dofs != length(prebasis)` here by design
@test length(get_prebasis(reffe)) == 18
@test Conformity(reffe) == DivConformity()
@test Pushforward(ArnoldWintherNC) == ReferenceFEs.DoubleContraVariantPiolaMap()
@test reffe == ReferenceFE(TRI, aw_nc, Float64)
# NB: `test_reference_fe` is not applicable here -- it asserts
# `num_dofs == length(prebasis)`, which an augmented element breaks by design.

# four DoFs per edge, three interior, none at the vertices
@test length.(get_face_own_dofs(reffe)) == [0, 0, 0, 4, 4, 4, 3]

@test evaluate(get_dof_basis(reffe), get_shapefuns(reffe)) ≈ Matrix(I, 15, 15)

@test_throws ErrorException ArnoldWintherNCRefFE(Float64, QUAD)
@test_throws ErrorException ArnoldWintherNCRefFE(Float64, TET)
@test_throws ErrorException ReferenceFE(TRI, aw_nc, Float64, 3)

############################################################################################
# The space is what the constraints say it is
############################################################################################

Ψ = get_shapefuns(reffe)

# (n⋅τn)|ₑ ∈ P₁(e): the degree-2 Legendre moment vanishes on every edge
μb2 = LegendreBasis(Val(1), Float64, 2)
quad_e = Quadrature(SEGMENT, 10)
ŝe = get_coordinates(quad_e)
ŵe = get_weights(quad_e)
for k in 1:3
  v1, v2 = get_face_coordinates(TRI, 1)[k]
  e = v2 - v1
  L = norm(e)
  n = ReferenceFEs._rot90(e / L)
  x = [v1 + s[1] * e for s in ŝe]
  Ψx = evaluate(Ψ, x)
  μ2 = evaluate(μb2, ŝe)[:, 3]
  for j in 1:15
    @test abs(sum(ŵe .* L .* [(n ⋅ ψ) ⋅ n for ψ in view(Ψx, :, j)] .* μ2)) < 1e-12
  end
end

# the space is a subspace of P₂(K;S) and contains P₁(K;S)
xs = get_coordinates(Quadrature(TRI, 6))
Ψxs = evaluate(Ψ, xs)
A = hcat([vcat([[ψ[1, 1], ψ[1, 2], ψ[2, 2]] for ψ in view(Ψxs, :, j)]...) for j in 1:15]...)
for τ in (x -> SymTensorValue{2,Float64}(1.0, 0.0, 0.0),
          x -> SymTensorValue{2,Float64}(0.0, 1.0, 0.0),
          x -> SymTensorValue{2,Float64}(0.0, 0.0, 1.0),
          x -> SymTensorValue{2,Float64}(x[1], 0.0, 0.0),
          x -> SymTensorValue{2,Float64}(0.0, x[2], 0.0),
          x -> SymTensorValue{2,Float64}(0.0, 0.0, x[1]))
  b = vcat([[τ(x)[1, 1], τ(x)[1, 2], τ(x)[2, 2]] for x in xs]...)
  @test norm(A * (A \ b) - b) < 1e-11
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
  P = copy(evaluate(FESpaces.AWNCChangeOfBasis(reffe, false), Jt, σ))
  Pinvt = copy(evaluate(FESpaces.AWNCChangeOfBasis(reffe, true), Jt, σ))

  @test transpose(Pinvt) * P ≈ Matrix(I, 15, 15)
  @test awnc_vandermonde(reffe, verts, σ) * P ≈ Matrix(I, 15, 15)

  if verts in (test_cells[2], test_cells[4])
    @test !isapprox(awnc_vandermonde(reffe, verts, σ), Matrix(I, 15, 15))
  end
end

############################################################################################
# FE space
############################################################################################

function test_awnc_fe_space(model)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 10)
  V = FESpace(model, ReferenceFE(TRI, aw_nc, Float64))

  topo = get_grid_topology(model)
  @test num_free_dofs(V) == 4*num_faces(topo, 1) + 3*num_cells(model)

  # P₁(K;S) ⊂ AWnc, so a global linear symmetric tensor is interpolated exactly
  τ(x) = SymTensorValue{2,Float64}(1.0 + 2*x[1] - x[2], 3.0 - x[1], 2.0 + x[2])
  τh = interpolate(τ, V)
  @test sqrt(sum(∫((τh - τ) ⊙ (τh - τ))dΩ)) < 1e-12

  w(x) = SymTensorValue{2,Float64}(sin(2*x[1]), cos(3*x[2]), sin(x[1] + x[2]))
  wh = interpolate(w, V)
  @test sqrt(sum(∫((wh - w) ⊙ (wh - w))dΩ)) > 1e-6

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 10)
  n = get_normal_vector(Λ).⁺

  # AWnc is nonconforming in H(div;S): the traction jumps ...
  @test sqrt(sum(∫((jump(wh) ⋅ n) ⋅ (jump(wh) ⋅ n))dΛ)) > 1e-4
  # ... but the normal-normal component is continuous outright, since its trace
  # is linear by construction and both of its moments are shared
  @test sqrt(sum(∫((n ⋅ jump(wh) ⋅ n) * (n ⋅ jump(wh) ⋅ n))dΛ)) < 1e-12
end

test_awnc_fe_space(unit_square(3))
test_awnc_fe_space(permute_cells(unit_square(3)))

############################################################################################
# The DIV operator
#
# `DIV` on a `DoubleContraVariantPiolaMap` space returns |det J| times the
# physical divergence -- the same normalization as the single Piola map -- so it
# pairs with a `ReferenceDomain()` measure.
############################################################################################

div_model = unit_square(3)
Ωd = Triangulation(div_model)
dΩd = Measure(Ωd, 10)
dωd = Measure(Ωd, 10, integration_domain_style=ReferenceDomain())

Vσ = FESpace(div_model, ReferenceFE(TRI, aw_nc, Float64))
Vu = FESpace(div_model, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1); conformity=:L2)

# a stress field inside the space, whose divergence is known in closed form
τd(x) = SymTensorValue{2,Float64}(1.0 + 2*x[1] - x[2], 3.0 - x[1], 2.0 + x[2])
divτd(x) = VectorValue(2.0, 0.0)
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
