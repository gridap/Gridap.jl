module MardalTaiWintherConvgTests

using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs, Gridap.Arrays, Gridap.TensorValues

# The Darcy-Stokes problem (u - ε²Δu + ∇p, div u) = (f, 0) of
# [Aznaran, Farrell & Kirby, SMAI-JCM 8 (2022) 399, (9.1)], for which the
# Mardal-Tai-Winther element is stable uniformly in ε, from the Stokes (ε ~ 1) to
# the Darcy (ε → 0) regime. The manufactured velocity is a curl, so it is
# divergence free and vanishes identically on ∂Ω.
#
# ∇(u) below is the *broken* gradient: the space is H¹-nonconforming.

function permute_cells(model, perms)
  grid = get_grid(model)
  cell_nodes = get_cell_node_ids(grid)
  permuted = [collect(cell_nodes[i])[perms[mod1(i, length(perms))]]
              for i in 1:length(cell_nodes)]
  pgrid = UnstructuredGrid(get_node_coordinates(grid), Table(permuted),
                           get_reffes(grid), get_cell_type(grid))
  pmodel = UnstructuredDiscreteModel(pgrid)
  topo = get_grid_topology(pmodel)
  labels = get_face_labeling(pmodel)
  D = num_cell_dims(pmodel)
  for d in 0:D-1
    get_face_entity(labels, d) .= ifelse.(Geometry.get_isboundary_face(topo, d), 2, 1)
  end
  get_face_entity(labels, D) .= 1
  add_tag!(labels, "interior", [1])
  add_tag!(labels, "boundary", [2])
  pmodel
end

permute_2d(m) = permute_cells(m, ([1,2,3],[2,3,1],[3,1,2],[2,1,3],[1,3,2],[3,2,1]))
permute_3d(m) = permute_cells(m, ([1,2,3,4],[2,3,4,1],[4,1,2,3],[2,1,3,4],[1,3,2,4],[4,3,2,1]))

function solve_darcy_stokes(model, ε, uex, pex, degree)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)
  f(x) = uex(x) - ε^2 * Δ(uex)(x) + ∇(pex)(x)

  V = FESpace(model, ReferenceFE(mtw, Float64, 1); dirichlet_tags="boundary")
  Q = FESpace(model, ReferenceFE(lagrangian, Float64, 0);
              conformity=:L2, constraint=:zeromean)
  X = MultiFieldFESpace([TrialFESpace(V, uex), TrialFESpace(Q)])
  Y = MultiFieldFESpace([V, Q])

  a((u, p), (v, q)) = ∫( u ⋅ v + ε^2 * (∇(u) ⊙ ∇(v)) - p * (∇ ⋅ v) - q * (∇ ⋅ u) )dΩ
  l((v, q)) = ∫( f ⋅ v )dΩ

  uh, ph = solve(AffineFEOperator(a, l, X, Y))
  (sqrt(sum(∫((uh - uex) ⋅ (uh - uex))dΩ)), sqrt(sum(∫((ph - pex) * (ph - pex))dΩ)))
end

function convg_test(models, hs, ε, uex, pex, degree)
  eu = Float64[]
  ep = Float64[]
  for model in models
    l2u, l2p = solve_darcy_stokes(model, ε, uex, pex, degree)
    println((l2u, l2p))
    push!(eu, l2u)
    push!(ep, l2p)
  end
  println(eu)
  println(ep)
  println(hs)
  eu, ep
end

function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end

############################################################################################
# 2D
############################################################################################

ψ(x) = (x[1]^2 * (1 - x[1])^2) * (x[2]^2 * (1 - x[2])^2)
u2(x) = VectorValue(∇(ψ)(x)[2], -∇(ψ)(x)[1])
p2(x) = cos(pi * x[1]) * cos(pi * x[2])

square(n) = simplexify(CartesianDiscreteModel((0,1,0,1), (n,n)))
ns2 = [4, 8, 16]
hs2 = [1/n for n in ns2]

println("--- 2D, second order in the velocity, first order in the P0 pressure ---")
for ε in (1.0, 1e-2, 1e-4)
  eu, ep = convg_test([square(n) for n in ns2], hs2, ε, u2, p2, 8)
  println("eps = $ε -- slope L2-norm u: $(slope(hs2,eu)), p: $(slope(hs2,ep))")
end
println("Sorted vs permuted: $(solve_darcy_stokes(square(8), 1e-4, u2, p2, 8)) vs " *
        "$(solve_darcy_stokes(permute_2d(square(8)), 1e-4, u2, p2, 8))")

############################################################################################
# 3D
#
# These meshes are strongly pre-asymptotic at ε = 1, so the rate is still
# climbing toward 2 there; in the Darcy limit the velocity is divergence free and
# the space is H(div)-conforming, so convergence is faster.
############################################################################################

Xt(t) = t^2 * (1 - t)^2
φ(x) = Xt(x[1]) * Xt(x[2]) * Xt(x[3])
u3(x) = VectorValue(∇(φ)(x)[2], -∇(φ)(x)[1], 0.0)
p3(x) = cos(pi * x[1]) * cos(pi * x[2]) * cos(pi * x[3])

cube(n) = simplexify(CartesianDiscreteModel((0,1,0,1,0,1), (n,n,n)))
ns3 = [2, 4, 8]
hs3 = [1/n for n in ns3]

println("--- 3D ---")
for ε in (1.0, 1e-4)
  eu, ep = convg_test([cube(n) for n in ns3], hs3, ε, u3, p3, 6)
  println("eps = $ε -- slope L2-norm u: $(slope(hs3,eu)), p: $(slope(hs3,ep))")
end
println("Sorted vs permuted: $(solve_darcy_stokes(cube(4), 1.0, u3, p3, 6)) vs " *
        "$(solve_darcy_stokes(permute_3d(cube(4)), 1.0, u3, p3, 6))")

end # module
