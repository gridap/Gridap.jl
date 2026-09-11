module HHJConvgTests

using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs, Gridap.Arrays

# The HHJ mixed method of [Arnold & Walker, SIAM J. Numer. Anal. 58 (2020) 2829,
# (2.4)]: find (σₕ, wₕ) ∈ Vₕ × Wₕ with
#
#   (σₕ, τ) + bₕ(τ, wₕ) = 0,   bₕ(σₕ, v) = -⟨f, v⟩,
#   bₕ(ϕ, v) = -Σ_T (ϕ, ∇²v)_T + Σ_E ⟨ϕ_nn, [∂v/∂n]⟩_E,
#
# Wₕ the continuous Lagrange space of degree r+1. The boundary edges are included
# in the second sum with [η] = η, which imposes ∂w/∂n = 0 weakly; together with
# w = 0 on ∂Ω, essential in Wₕ, this is the clamped plate. The manufactured
# solution is a product of x²(1-x)², which vanishes together with its gradient
# on ∂Ω.

function solve_plate(nc, r, wex, σex, f; permute=false)
  model = simplexify(CartesianDiscreteModel((0,1,0,1), nc))
  permute && (model = permute_cells(model))
  deg = 2*r + 6
  Ω = Triangulation(model)
  dΩ = Measure(Ω, deg)
  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, deg)
  nΛ = get_normal_vector(Λ).⁺
  Γ = BoundaryTriangulation(model)
  dΓ = Measure(Γ, deg)
  nΓ = get_normal_vector(Γ)

  Vσ = FESpace(model, ReferenceFE(hhj, Float64, r))
  Vw = FESpace(model, ReferenceFE(lagrangian, Float64, r + 1); dirichlet_tags="boundary")
  X = MultiFieldFESpace([Vσ, TrialFESpace(Vw, 0.0)])
  Y = MultiFieldFESpace([Vσ, Vw])

  bh(ϕ, v) = ∫( -(ϕ ⊙ ∇∇(v)) )dΩ +
             ∫( (nΛ ⋅ mean(ϕ) ⋅ nΛ) * (jump(∇(v)) ⋅ nΛ) )dΛ +
             ∫( (nΓ ⋅ ϕ ⋅ nΓ) * (∇(v) ⋅ nΓ) )dΓ
  A((σ, w), (τ, v)) = ∫( σ ⊙ τ )dΩ + bh(τ, w) + bh(σ, v)
  L((τ, v)) = ∫( -(f * v) )dΩ

  σh, wh = solve(AffineFEOperator(A, L, X, Y))

  ew = wh - wex
  eσ = σh - σex
  l2w = sqrt( sum( ∫( ew * ew )dΩ ) )
  l2σ = sqrt( sum( ∫( eσ ⊙ eσ )dΩ ) )

  return l2w, l2σ
end

# The same mesh with the vertices of each cell listed in a permuted order: the
# discrete space must not depend on it. For r ≥ 1 this exercises the sign flip of
# the odd-degree edge DoFs.
function permute_cells(model, perms=([1,2,3],[2,3,1],[3,1,2],[2,1,3],[1,3,2],[3,2,1]))
  grid = get_grid(model)
  cell_nodes = get_cell_node_ids(grid)
  permuted = [collect(cell_nodes[i])[perms[mod1(i, length(perms))]]
              for i in 1:length(cell_nodes)]
  pgrid = UnstructuredGrid(get_node_coordinates(grid), Table(permuted),
                           get_reffes(grid), get_cell_type(grid))
  pmodel = UnstructuredDiscreteModel(pgrid)
  topo = get_grid_topology(pmodel)
  labels = get_face_labeling(pmodel)
  for d in 0:1
    get_face_entity(labels, d) .= ifelse.(Geometry.get_isboundary_face(topo, d), 2, 1)
  end
  get_face_entity(labels, 2) .= 1
  add_tag!(labels, "interior", [1])
  add_tag!(labels, "boundary", [2])
  pmodel
end

function convg_test(ncs, r, wex, σex, f)
  ew = Float64[]
  eσ = Float64[]
  hs = Float64[]
  for nc in ncs
    l2w, l2σ = solve_plate(nc, r, wex, σex, f)
    println((l2w, l2σ))
    push!(ew, l2w)
    push!(eσ, l2σ)
    h = 1/nc[1]
    push!(hs, h)
  end
  println(ew)
  println(eσ)
  println(hs)
  ew, eσ, hs
end

function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end

p(x) = x^2 * (1 - x)^2
pdd(x) = 2 - 12*x + 12*x^2
w(x) = p(x[1]) * p(x[2])
σ = ∇∇(w)
f(x) = 24*p(x[2]) + 2*pdd(x[1]) * pdd(x[2]) + 24*p(x[1])

# O(h^{r+2}) in the displacement, O(h^{r+1}) in the moment tensor
ncs = [(4,4),(8,8),(16,16)]

for r in 0:1
  ew, eσ, hs = convg_test(ncs, r, w, σ, f)
  println("Order $r -- slope L2-norm w: $(slope(hs,ew))")
  println("Order $r -- slope L2-norm σ: $(slope(hs,eσ))")

  sorted = solve_plate((8,8), r, w, σ, f)
  permuted = solve_plate((8,8), r, w, σ, f; permute=true)
  println("Order $r -- sorted vs permuted mesh: $(sorted) vs $(permuted)")
end

end # module
