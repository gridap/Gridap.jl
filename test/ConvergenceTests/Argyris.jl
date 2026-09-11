module ArgyrisConvgTests

using Gridap
using Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs, Gridap.Arrays

# The clamped plate Δ²u = f, discretized with the conforming Argyris form
# a(u,v) = ∫ D²u : D²v -- no broken norms, the space is C¹. The boundary DoFs are
# set by interpolating the exact solution: it does not satisfy D²u = 0 there, and
# constraining the boundary Hessian to its exact values is consistent.

function solve_biharmonic(nc, uex, Huex, f; permute=false)
  model = simplexify(CartesianDiscreteModel((0,1,0,1), nc))
  permute && (model = permute_cells(model))
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 12)

  V = FESpace(model, ReferenceFE(argyris, Float64, 5); dirichlet_tags="boundary")
  U = TrialFESpace(V, uex)

  a(u, v) = ∫( ∇∇(u) ⊙ ∇∇(v) )dΩ
  l(v) = ∫( f * v )dΩ

  uh = solve(AffineFEOperator(a, l, U, V))

  eu = uh - uex
  eH = ∇∇(uh) - Huex
  l2u = sqrt( sum( ∫( eu * eu )dΩ ) )
  h2u = sqrt( sum( ∫( eH ⊙ eH )dΩ ) )

  return l2u, h2u
end

# The same mesh with the vertices of each cell listed in a permuted order: the
# discrete space must not depend on it. Both edge orientations occur, since the
# list holds even and odd permutations.
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

function convg_test(ncs, uex, Huex, f)
  el2 = Float64[]
  eh2 = Float64[]
  hs = Float64[]
  for nc in ncs
    l2, h2 = solve_biharmonic(nc, uex, Huex, f)
    println((l2, h2))
    push!(el2, l2)
    push!(eh2, h2)
    h = 1/nc[1]
    push!(hs, h)
  end
  println(el2)
  println(eh2)
  println(hs)
  el2, eh2, hs
end

function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end

p(x) = x^2 * (1 - x)^2
pdd(x) = 2 - 12*x + 12*x^2
u(x) = p(x[1]) * p(x[2])
Hu = ∇∇(u)
f(x) = 24*p(x[2]) + 2*pdd(x[1]) * pdd(x[2]) + 24*p(x[1])

# Not pushed further: the L2 error is already ~3e-8 at 8x8, and the Argyris
# system is badly conditioned without the DoF rescaling of
# [Aznaran, Farrell & Kirby, SMAI-JCM 8 (2022) 399, §5.4], which is not implemented.
ncs = [(2,2),(4,4),(8,8)]

el, eh, hs = convg_test(ncs, u, Hu, f)
println("Slope L2-norm u: $(slope(hs,el))")
println("Slope H2-seminorm u: $(slope(hs,eh))")

# The solved solution must not depend on how the cells list their vertices. The
# boundary data is interpolated, so this also exercises the discretized edge
# moment; here the internal quadrature integrates it exactly, so the agreement is
# to round-off.
sorted = solve_biharmonic((4,4), u, Hu, f)
permuted = solve_biharmonic((4,4), u, Hu, f; permute=true)
println("Sorted vs permuted mesh: $(sorted) vs $(permuted)")

end # module
