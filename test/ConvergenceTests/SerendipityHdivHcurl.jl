module SerendipityHdivHcurlConvgTests

using Test
using Random
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs

# The same mesh with the vertices of each cell relabelled by a symmetry of the
# cube: the discrete solution must not depend on it
function relabel_cells(model)
  grid = get_grid(model)
  syms = get_vertex_permutations(only(get_polytopes(model)))
  perms = syms[randperm(MersenneTwister(0), length(syms))]
  cell_nodes = get_cell_node_ids(grid)
  relabelled = [collect(cell_nodes[i])[perms[mod1(i,length(perms))]] for i in 1:length(cell_nodes)]
  x = vec(collect(get_node_coordinates(grid)))
  rmodel = UnstructuredDiscreteModel(UnstructuredGrid(
    x, Table(relabelled), collect(get_reffes(grid)), collect(get_cell_type(grid))))
  topo = get_grid_topology(rmodel)
  labels = get_face_labeling(rmodel)
  D = num_cell_dims(rmodel)
  for d in 0:D-1
    get_face_entity(labels, d) .= ifelse.(Geometry.get_isboundary_face(topo, d), 2, 1)
  end
  get_face_entity(labels, D) .= 1
  add_tag!(labels, "interior", [1])
  add_tag!(labels, "boundary", [2])
  rmodel
end

l2(u, dΩ) = sqrt(sum(∫(u⋅u)dΩ))
rates(errs) = [log2(errs[i]/errs[i+1]) for i in 1:length(errs)-1]

# Mixed Poisson, u = -∇p, div u = f, p = 0 on ∂Ω, with the broken 𝓟_q pressure
function mixed_poisson_errors(model, reffe, q)
  p(x) = prod(sin(π*xi) for xi in Tuple(x))
  u(x) = -∇(p)(x)
  f(x) = (∇⋅u)(x)
  V = FESpace(model, reffe)
  Q = FESpace(model, ReferenceFE(lagrangian, Float64, q; space=:P), conformity=:L2)
  X = MultiFieldFESpace([V, Q])
  dΩ = Measure(Triangulation(model), 2*get_order(reffe)+4)
  a((uh, ph), (v, w)) = ∫(uh⋅v - ph*(∇⋅v) + (∇⋅uh)*w)dΩ
  l((v, w)) = ∫(f*w)dΩ
  uh, ph = solve(AffineFEOperator(a, l, X, X))
  (l2(u - uh, dΩ), l2(p - ph, dΩ))
end

# curl curl u + u = f, u×n prescribed on ∂Ω
function maxwell_errors(model, reffe)
  D = num_cell_dims(model)
  u(x) = D == 2 ? VectorValue(sin(π*x[2])*exp(x[1]), cos(π*x[1])*x[2]^2) :
                  VectorValue(sin(π*x[2])*exp(x[3]), cos(π*x[3])*x[1]^2, sin(π*x[1])*x[2]^2)
  rot(s) = x -> (g = ∇(s)(x); VectorValue(g[2], -g[1])) # curl of a scalar, in 2D
  f(x) = (D == 2 ? rot(∇×u)(x) : (∇×(∇×u))(x)) + u(x)
  V = FESpace(model, reffe, dirichlet_tags="boundary")
  U = TrialFESpace(V, u)
  dΩ = Measure(Triangulation(model), 2*get_order(reffe)+4)
  a(uh, v) = ∫((∇×uh)⋅(∇×v) + uh⋅v)dΩ
  l(v) = ∫(f⋅v)dΩ
  uh = solve(AffineFEOperator(a, l, U, V))
  (l2(u - uh, dΩ), l2((∇×u) - (∇×uh), dΩ))
end

square(n) = CartesianDiscreteModel((0,1,0,1), (n,n))
cube(n) = CartesianDiscreteModel((0,1,0,1,0,1), (n,n,n))

# 𝓢ᵣΛ¹ × 𝓟ᵣ₋₁: ‖u - uₕ‖ = O(h^{r+1}), ‖p - pₕ‖ = O(hʳ)
# 𝓢⁻ₖ₊₁Λ¹ × 𝓟ₖ: both O(h^{k+1})
for (name, o, q, ru, rp) in ((sdiv,1,0,2,1), (sdiv,2,1,3,2), (trimmed_sdiv,0,0,1,1), (trimmed_sdiv,1,1,2,2))
  reffe = ReferenceFE(QUAD, name, Float64, o)
  errs = [mixed_poisson_errors(relabel_cells(square(n)), reffe, q) for n in (4,8,16)]
  println("$name, order $o: ", errs)
  @test all(rates(first.(errs)) .> ru - 0.2)
  @test all(rates(last.(errs)) .> rp - 0.2)
  @test collect(mixed_poisson_errors(square(8), reffe, q)) ≈ collect(errs[2]) rtol=1e-8
end

# 𝓢ᵣΛ¹: ‖u - uₕ‖ = O(h^{r+1}), ‖curl(u - uₕ)‖ = O(hʳ); 𝓢⁻ₖ₊₁Λ¹: both O(h^{k+1})
for (name, o, ru, rc) in ((scurl,1,2,1), (scurl,2,3,2), (trimmed_scurl,0,1,1), (trimmed_scurl,1,2,2))
  reffe = ReferenceFE(QUAD, name, Float64, o)
  errs = [maxwell_errors(relabel_cells(square(n)), reffe) for n in (4,8,16)]
  println("$name, order $o: ", errs)
  @test all(rates(first.(errs)) .> ru - 0.2)
  @test all(rates(last.(errs)) .> rc - 0.2)
  @test collect(maxwell_errors(square(8), reffe)) ≈ collect(errs[2]) rtol=1e-8
end

# 3D, on sorted and relabelled meshes
for name in (sdiv, trimmed_sdiv)
  reffe = ReferenceFE(HEX, name, Float64, 1)
  @test collect(mixed_poisson_errors(cube(2), reffe, name == sdiv ? 0 : 1)) ≈
        collect(mixed_poisson_errors(relabel_cells(cube(2)), reffe, name == sdiv ? 0 : 1)) rtol=1e-8
end
for name in (scurl, trimmed_scurl)
  reffe = ReferenceFE(HEX, name, Float64, 1)
  @test collect(maxwell_errors(cube(2), reffe)) ≈ collect(maxwell_errors(relabel_cells(cube(2)), reffe)) rtol=1e-8
end

end # module
