module ModelWithFEMapsTests

using Test
using Gridap
# using Gridap.Arrays
# using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
# using Gridap.Geometry: DiscreteModelMock
# using Gridap.Io
# using Gridap.CellData

# using Gridap.FESpaces

# include("ModelWithFEMaps.jl")

# Basic tests

domain = (0,1)
partition = (2,)
model = CartesianDiscreteModel(domain,partition)

T = VectorValue{1,Float64}; order = 1; pol = Polytope(HEX_AXIS...)
reffe = LagrangianRefFE(T,pol,order)

Vₕ = FESpace(model,reffe;conformity=:H1)
u(x) = 2*x
Uₕ = TrialFESpace(Vₕ,u)
uh = interpolate_everywhere(u,Uₕ)

geo_map = uh

new_model = MappedDiscreteModel(model,uh)

grid = get_grid(model)
new_grid = MappedGrid(grid,uh)
new_grid_2 = get_grid(new_model)

get_node_coordinates(new_grid.grid)
get_node_coordinates(new_model)
new_grid.geo_map

parent_node_to_coords = get_node_coordinates(new_grid.grid)
lazy_map(new_grid.geo_map,parent_node_to_coords)

get_node_coordinates(grid)
collect(get_node_coordinates(new_grid))
get_node_coordinates(new_grid_2)

@test get_node_coordinates(new_grid) == get_node_coordinates(grid) * 2
@test get_node_coordinates(new_grid) == get_node_coordinates(new_grid_2)
@test get_node_coordinates(new_grid) == get_node_coordinates(new_model)

Tₕ = Triangulation(new_model)
@test get_node_coordinates(Tₕ) == get_node_coordinates(new_grid)

test_discrete_model(new_model)
test_grid(grid)
test_triangulation(Tₕ)

# FE problem test


domain = (0,1)
partition = (2)
model = CartesianDiscreteModel(domain,partition)

function get_map_model()

  T = VectorValue{1,Float64}; order = 1; pol = Polytope(HEX_AXIS...)
  reffe = LagrangianRefFE(T,pol,order)

  Vₕ = FESpace(model,reffe;conformity=:H1)
  d(x) = 2*x            # Analytical solution (for Dirichlet data)
  Uₕ = TrialFESpace(Vₕ,d)
  dh = interpolate_everywhere(d,Uₕ)
  new_model = MappedModel(model,dh)
end

map_model = get_map_model()

T = Float64; order = 1; pol = Polytope(HEX_AXIS...)
reffe = LagrangianRefFE(T,pol,order)

Vₕ = FESpace(map_model,reffe;conformity=:H1,dirichlet_tags="boundary")
u(x) = 2*x[1]            # Analytical solution (for Dirichlet data)
Uₕ = TrialFESpace(Vₕ,u)
uh = interpolate_everywhere(u,Uₕ)
get_free_dof_values(uh)

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
@test get_free_dof_values(uh_map) == 2*get_free_dof_values(uh)

domain = (0,2)
partition = (2)
model = CartesianDiscreteModel(domain,partition)
uh = fe_problem(model)
@test get_free_dof_values(uh_map) == get_free_dof_values(uh)

end # module
