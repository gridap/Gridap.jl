module HermiteConvgTests

using Gridap
using Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs, Gridap.Arrays

# The Poisson problem -Δu = f with the cubic Hermite element, the test problem
# of [Kirby, SMAI-JCM 4 (2018) 197, §5.3]: fourth order in L² and third in H¹,
# as for cubic Lagrange. The boundary DoFs are set by interpolating the exact
# solution, which also fixes the boundary gradient DoFs -- the tangential ones
# are what a Dirichlet condition implies, the normal ones are a consistent
# extra constraint.

function solve_poisson(nc, uex, f; permute=false)
  model = simplexify(CartesianDiscreteModel((0,1,0,1), nc))
  permute && (model = permute_cells(model))
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 8)

  V = FESpace(model, ReferenceFE(hermite, Float64, 3); dirichlet_tags="boundary")
  U = TrialFESpace(V, uex)

  a(u, v) = ∫( ∇(u) ⋅ ∇(v) )dΩ
  l(v) = ∫( f * v )dΩ

  uh = solve(AffineFEOperator(a, l, U, V))

  eu = uh - uex
  eg = ∇(uh) - ∇(uex)
  l2u = sqrt( sum( ∫( eu * eu )dΩ ) )
  h1u = sqrt( sum( ∫( eg ⋅ eg )dΩ ) )

  return l2u, h1u
end

# The same mesh with the vertices of each cell listed in a permuted order: the
# discrete space must not depend on it.
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

function convg_test(ncs, uex, f)
  el2 = Float64[]
  eh1 = Float64[]
  hs = Float64[]
  for nc in ncs
    l2, h1 = solve_poisson(nc, uex, f)
    println((l2, h1))
    push!(el2, l2)
    push!(eh1, h1)
    h = 1/nc[1]
    push!(hs, h)
  end
  println(el2)
  println(eh1)
  println(hs)
  el2, eh1, hs
end

function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end

u(x) = sin(2π*x[1]) * sin(2π*x[2])
f(x) = 8π^2 * u(x)

ncs = [(4,4),(8,8),(16,16),(32,32)]

el, eh, hs = convg_test(ncs, u, f)
println("Slope L2-norm u: $(slope(hs,el))")
println("Slope H1-seminorm u: $(slope(hs,eh))")

# The solved solution must not depend on how the cells list their vertices.
sorted = solve_poisson((8,8), u, f)
permuted = solve_poisson((8,8), u, f; permute=true)
println("Sorted vs permuted mesh: $(sorted) vs $(permuted)")

end # module
