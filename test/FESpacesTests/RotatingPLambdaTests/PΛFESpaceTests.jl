module PΛFESpaceTests
# Mesh-level conformity layer: H(curl)-conforming FESpaces on simplicial
# meshes with ARBITRARY local vertex orderings, for BOTH the rotating P_rΛ¹
# and the trimmed P_r⁻Λ¹ bases. The per-cell change of basis is the rotation
# calculus at π_K = sortperm(cell global vertex ids); no per-face data, no
# sign flips, identity face-own-dof permutations for every pindex.
#
# These tests PIN the convention implemented by compute_pλ_change
# (cell_change[i,j] = C(π_K)[s(j),i], cell_change_invt[i,j] =
# C(invperm(π_K))[i,s(j)], s the face-wise sorted relabelling): during
# development the competing candidates (transposes, inverse permutation,
# un-relabelled variants) were iterated programmatically and each fails the
# machine-zero tangential-jump assertions below at O(1).

using Gridap
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Fields
using Gridap.Arrays: Table, Fill
using Combinatorics: permutations
using LinearAlgebra
using LinearAlgebra: eigvals, issymmetric
using Test

# ─── Mesh builders ───────────────────────────────────────────────────────────
# Two D-simplices sharing the facet with global vertices {2,…,D+1}; the cells'
# local vertex orderings are arbitrary permutations of their vertex sets.

const X2D = [[0.0,0.0], [1.0,0.0], [0.0,1.0], [1.0,1.0]]
const X3D = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0], [1.0,1.0,1.0]]

function two_cell_model(D::Int, c1::Vector{Int}, c2::Vector{Int})
  X = D == 2 ? X2D : X3D
  coords = [Point(x...) for x in X]
  reffe  = LagrangianRefFE(Float64, D == 2 ? TRI : TET, 1)
  grid   = UnstructuredGrid(coords, Table([c1, c2]), [reffe], fill(Int8(1), 2))
  UnstructuredDiscreteModel(grid)
end

# Affine geometry (o, J) of cell c: F(x̂) = o + J x̂.
function cell_affine(model, c)
  X    = get_node_coordinates(get_grid(model))
  conn = get_cell_node_ids(get_grid(model))
  g = conn[c]
  D = length(first(X))
  o = collect(X[g[1]].data)
  J = hcat([collect(X[g[j+1]].data) .- o for j in 1:D]...)
  o, J
end

# Physical points on the shared facet {2,…,D+1} and its tangent vectors.
function interface_points_tangents(D::Int)
  X = D == 2 ? X2D : X3D
  ts = [X[2+j] .- X[2] for j in 1:D-1]
  wts = D == 2 ? [(0.15,), (0.5,), (0.83,)] :
                 [(0.2, 0.3), (0.5, 0.25), (0.1, 0.7), (1/3, 1/3)]
  pts = [X[2] .+ sum(w[j] .* ts[j] for j in 1:D-1) for w in wts]
  pts, ts
end

# ─── Trace jump across the shared facet ──────────────────────────────────────
# cell_fields[c] must be a vector of fields on the reference cell returning
# physical covector values (the FESpace's own cell basis, or a manually pushed
# one). Pairs local dofs via cell_ids; returns the max over interface points
# and tangent directions of
#   • |trace from cell 1 − trace from cell 2| for every shared global dof,
#   • |trace| for every non-shared dof (must vanish on the facet).
function interface_jump(model, cell_fields, cell_ids, D)
  pts, ts = interface_points_tangents(D)
  d = [collect(cell_ids[1]), collect(cell_ids[2])]
  shared = intersect(Set(d[1]), Set(d[2]))
  geo = [cell_affine(model, c) for c in 1:2]
  maxjump = 0.0
  for p in pts
    tr = Dict{Int,Vector{Vector{Float64}}}()
    for c in 1:2
      o, J = geo[c]
      x̂ = Point((J \ (p .- o))...)
      vals = evaluate(cell_fields[c], [x̂])
      for (l, gdof) in enumerate(d[c])
        v = collect(vals[1, l].data)
        tv = [dot(v, t) for t in ts]
        if gdof in shared
          push!(get!(tr, gdof, Vector{Float64}[]), tv)
        else
          maxjump = max(maxjump, maximum(abs, tv))
        end
      end
    end
    for (_, vs) in tr
      maxjump = max(maxjump, maximum(abs, vs[1] .- vs[2]))
    end
  end
  maxjump
end

fespace_jump(model, V, D) =
  interface_jump(model, Gridap.CellData.get_data(get_fe_basis(V)),
                 get_cell_dof_ids(V), D)

# ─── Pointwise interpolation round trip ──────────────────────────────────────
function interpolation_error(model, V, u, D)
  uh  = interpolate(u, V)
  uhd = Gridap.CellData.get_data(uh)
  x̂s  = D == 2 ? [Point(0.2, 0.3), Point(0.5, 0.25), Point(0.1, 0.7)] :
                 [Point(0.1, 0.2, 0.3), Point(0.25, 0.25, 0.25), Point(0.05, 0.5, 0.1)]
  err = 0.0
  for c in 1:2
    o, J = cell_affine(model, c)
    for x̂ in x̂s
      p = Point((o .+ J * collect(x̂.data))...)
      v = evaluate(uhd[c], x̂)
      err = max(err, maximum(abs, collect(v.data) .- collect(u(p).data)))
    end
  end
  err
end

# Interpolable fields per space
# constants: all r (full & trimmed)
u_const(D) = D == 2 ? (x -> VectorValue(0.7, -1.3)) :
                      (x -> VectorValue(0.7, -1.3, 0.4))
# affine:  any r (full),  r ≥ 2 (trimmed)
u_lin(D) = D == 2 ?
  (x -> VectorValue(1.0 + 2x[1] - x[2], -0.5 + x[1] + 3x[2])) :
  (x -> VectorValue(1.0 + 2x[1] - x[2] + x[3],
                    -0.5 + x[1] + 3x[2] - 2x[3],
                    0.25 - x[1] + x[2] + x[3]))
# (Whitney) the Koszul part κ(dx¹∧dx²) = x¹dx² − x²dx¹  for trimmed r = 1
u_koszul(D) = D == 2 ? (x -> VectorValue(-x[2], x[1])) :
                       (x -> VectorValue(-x[2], x[1], 0.0))

# Exterior calculus form of the same fields
#
# u_const(D) = D == 2 ? (x -> DifferentialFormValue{1,2}((0.7, -1.3))) :
#                       (x -> DifferentialFormValue{1,3}((0.7, -1.3, 0.4)))
# u_lin(D) = D == 2 ?
#   (x -> DifferentialFormValue{1,2}((1.0 + 2x[1] - x[2], -0.5 + x[1] + 3x[2]))) :
#   (x -> DifferentialFormValue{1,3}((1.0 + 2x[1] - x[2] + x[3],
#                                     -0.5 + x[1] + 3x[2] - 2x[3],
#                                     0.25 - x[1] + x[2] + x[3])))
# u_koszul(D) = D == 2 ? (x -> DifferentialFormValue{1,2}((-x[2], x[1]))) :
#                        (x -> DifferentialFormValue{1,3}((-x[2], x[1], 0.0)))

test_fields(name, r, D) =
  name == :trimmed && r == 1 ? [u_const(D), u_koszul(D)] : [u_const(D), u_lin(D)]

simplex(D) = D == 2 ? TRI : TET
make_reffe(name, D, r) = name == :rotating ? RotatingPΛRefFE(Float64, simplex(D), r) :
                                             TrimmedPΛRefFE(Float64, simplex(D), r)

# ─── Two-triangle conformity: all relative orderings, both bases, r = 1,2,3 ──

const TRI1_ORDERINGS = ([1,2,3], [3,1,2], [2,1,3])

@testset "2D conformity + interpolation, all orderings" begin
  for name in (:rotating, :trimmed), r in (1, 2, 3)
    rf = make_reffe(name, 2, r)
    us = test_fields(name, r, 2)
    @testset "$name r=$r" begin
      for c1 in TRI1_ORDERINGS, c2 in permutations([2, 3, 4])
        model = two_cell_model(2, c1, collect(c2))
        V = FESpace(model, rf, conformity=:HCurl)
        @test fespace_jump(model, V, 2) < 1e-12
        for u in us
          @test interpolation_error(model, V, u, 2) < 1e-10
        end
      end
    end
  end
end

# ─── Two-tet conformity: scrambles incl. non-cyclic shared-face perm ─────────
# Shared face = global vertices {2,3,4}. [3,2,4,5] lists it with the
# NON-CYCLIC (odd) permutation (3,2,4) of (2,3,4).

const TET_ORDERINGS = (
  ([1,2,3,4], [2,3,4,5]),   # both sorted (fast path, no change of basis)
  ([1,2,3,4], [3,2,4,5]),   # non-cyclic permutation of the shared face
  ([1,2,3,4], [4,3,2,5]),
  ([1,2,3,4], [5,4,3,2]),
  ([2,4,1,3], [5,4,3,2]),
  ([4,1,3,2], [3,2,5,4]),
)

@testset "3D conformity + interpolation, scrambled tets" begin
  for name in (:rotating, :trimmed), r in (1, 2)
    rf = make_reffe(name, 3, r)
    us = test_fields(name, r, 3)
    @testset "$name r=$r" begin
      for (c1, c2) in TET_ORDERINGS
        model = two_cell_model(3, c1, c2)
        V = FESpace(model, rf, conformity=:HCurl)
        @test fespace_jump(model, V, 3) < 1e-12
        for u in us
          @test interpolation_error(model, V, u, 3) < 1e-10
        end
      end
    end
  end
end

# ─── Fast path: sorted cells produce no change of basis ──────────────────────

@testset "sorted mesh takes the nothing fast path" begin
  model = two_cell_model(2, [1,2,3], [2,3,4])
  rf = RotatingPΛRefFE(Float64, TRI, 1)
  cell_reffe = Fill(rf, 2)
  name = ReferenceFEs.get_name(eltype(cell_reffe))
  ch = FESpaces.compute_cell_bases_changes(
    name, CoVariantPiolaMap(), model, cell_reffe, nothing)
  @test ch === nothing
  # and the resulting space is still conforming
  V = FESpace(model, rf, conformity=:HCurl)
  @test fespace_jump(model, V, 2) < 1e-12
end

# ─── Negative control: bypassing the rotation breaks conformity ──────────────
# Same mesh, same dof gluing, but raw covariant Piola pushforward with
# cell_changes = nothing: the tangential traces must NOT match (O(1) jump).

@testset "negative control: no rotation → O(1) jump" begin
  for name in (:rotating, :trimmed)
    rf = make_reffe(name, 2, 2)
    model = two_cell_model(2, [1,2,3], [4,3,2])
    V = FESpace(model, rf, conformity=:HCurl)          # for the dof pairing
    @test fespace_jump(model, V, 2) < 1e-12            # rotation path is fine
    cell_reffe = Fill(rf, 2)
    cell_Jt = lazy_map(Broadcasting(∇), get_cell_map(get_grid(model)))
    cell_fields_raw, _ = FESpaces.get_cell_shapefuns_and_dof_basis(
      CoVariantPiolaMap(), model, cell_reffe, nothing, cell_Jt)
    @test interface_jump(model, cell_fields_raw, get_cell_dof_ids(V), 2) > 0.1
  end
end

# ─── Assembly smoke tests: SPD mass matrix + curl-curl (gradient path) ───────

function assembly_checks(model, rf; deg=6)
  trian = Triangulation(model)
  dΩ    = Measure(trian, deg)
  V     = FESpace(model, rf, conformity=:HCurl)
  U     = TrialFESpace(V)
  n     = num_free_dofs(V)
  # ⟨u,v⟩ and ⟨du,dv⟩ of the 1-forms, as the vector proxy identities
  # vol_coeff(u ∧ ⋆v) = u⋅v and vol_coeff(du ∧ ⋆dv) = (∇×u)⋅(∇×v).
  a_mass(u, v)  = ∫(u ⋅ v) * dΩ
  a_stiff(u, v) = ∫((∇×u) ⋅ (∇×v)) * dΩ
  # a_mass(u, v)  = ∫(vol_coeff(u ∧ hodge_star(v))) * dΩ
  # a_stiff(u, v) = ∫(vol_coeff(d_1form(u) ∧ hodge_star(d_1form(v)))) * dΩ
  M = assemble_matrix(a_mass, U, V)
  A = assemble_matrix(a_stiff, U, V)
  n, Matrix(M), Matrix(A)
end

@testset "assembly on scrambled meshes: mass SPD, curl-curl sym PSD" begin
  cases = (
    (:rotating, 2, 2, ([3,1,2], [4,3,2])),
    (:trimmed,  2, 2, ([3,1,2], [4,3,2])),
    (:rotating, 3, 1, ([2,4,1,3], [5,4,3,2])),
    (:trimmed,  3, 2, ([1,2,3,4], [3,2,4,5])),
  )
  for (name, D, r, (c1, c2)) in cases
    model = two_cell_model(D, c1, c2)
    n, M, A = assembly_checks(model, make_reffe(name, D, r))
    @test n > 0
    @test size(M) == (n, n)
    @test issymmetric(M)
    @test isposdef(M)
    @test issymmetric(A)
    @test all(e -> e > -1e-10, real.(eigvals(A)))   # curl has a null space
  end
end

@testset "assembly on a simplexified Cartesian mesh (sorted fast path)" begin
  model = simplexify(CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (2, 2)))
  for name in (:rotating, :trimmed)
    trian = Triangulation(model)
    dΩ    = Measure(trian, 4)
    V     = FESpace(model, make_reffe(name, 2, 1), conformity=:HCurl)
    U     = TrialFESpace(V)
    a_mass(u, v) = ∫(u ⋅ v) * dΩ
    # a_mass(u, v) = ∫(vol_coeff(u ∧ hodge_star(v))) * dΩ
    M = Matrix(assemble_matrix(a_mass, U, V))
    @test issymmetric(M)
    @test isposdef(M)
  end
end

end # module
