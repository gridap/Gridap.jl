module MappedDiscreteModelsTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData


# Basic tests
domain = (0,1)
partition = (2,)
model = CartesianDiscreteModel(domain,partition)

u(x) = 2*x

grid = get_grid(model)
new_grid = MappedGrid(grid,u)

nodes = get_node_coordinates(new_grid.grid)
new_nodes = get_node_coordinates(new_grid)

@test nodes*2 ≈ new_nodes
@test u.(nodes) ≈ new_nodes



new_model = MappedDiscreteModel(model,u)
new_grid = MappedGrid(get_grid(model),u)
new_grid_2 = get_grid(new_model)


@test get_node_coordinates(new_grid) ≈ get_node_coordinates(grid) * 2
@test get_node_coordinates(new_grid) ≈ get_node_coordinates(new_grid_2)
@test get_node_coordinates(new_grid) ≈ get_node_coordinates(new_model)

Tₕ = Triangulation(new_model)
@test get_node_coordinates(Tₕ) == get_node_coordinates(new_grid)

test_discrete_model(new_model)
test_grid(new_grid)
test_triangulation(Tₕ)

# FE problem test

u(x) = 2*x
domain = (0,1)
partition = (2)
model = CartesianDiscreteModel(domain,partition)
new_model = MappedDiscreteModel(model,u)

function fe_problem(model)

  T = Float64; order = 1; pol = Polytope(HEX_AXIS...)
  reffe = LagrangianRefFE(T,pol,order)

  Vₕ = FESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
  u(x) = 2*x[1]            # Analytical solution (for Dirichlet data)
  Uₕ = TrialFESpace(Vₕ,u)
  uh = interpolate_everywhere(u,Uₕ)


  degree = 2
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  f = 0.0
  a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ
  b(v) = ∫( v*f )*dΩ

  op = AffineFEOperator(a,b,Uₕ,Vₕ)

  return uh = solve(op)
end

uh_map = fe_problem(new_model)
uh = fe_problem(model)
@test get_free_dof_values(uh_map) ≈ 2*get_free_dof_values(uh)

domain = (0,2)
partition = (2)
model = CartesianDiscreteModel(domain,partition)
uh = fe_problem(model)
@test get_free_dof_values(uh_map) ≈ get_free_dof_values(uh)

end # module
