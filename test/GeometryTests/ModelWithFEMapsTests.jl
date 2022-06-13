module ModelWithFEMapsTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Geometry: DiscreteModelMock
using Gridap.Io
using Gridap.CellData

using Gridap.FESpaces

struct GridWithFEMap{Dc,Dp,T} <: Gridap.Geometry.Grid{Dc,Dp}
  grid::Grid{Dc,Dp}
  geo_map::T
  function GridWithFEMap(grid::Grid{Dc,Dp},geo_map) where {Dc,Dp}
    model_map=get_cell_map(grid)
    map=lazy_map(∘,geo_map,model_map)
    new{Dc,Dp,typeof(map)}(model,map)
  end
end

function GridWithFEMap(grid::Grid{Dc,Dp},uh::FEFunction) where {Dc,Dp}
  GridWithFEMap(grid,Gridap.CellData.get_data(uh))
end

import Gridap.ReferenceFEs.get_node_coordinates
import Gridap.Geometry.get_cell_node_ids
import Gridap.Geometry.get_reffes
import Gridap.Geometry.get_cell_type
import Gridap.Geometry.Grid

function get_node_coordinates(grid::GridWithFEMap)
  parent_node_to_coords = get_node_coordinates(grid.grid)
  lazy_map(geo_map,parent_node_to_coords)
end

get_cell_node_ids(grid::GridWithFEMap) = get_cell_node_ids(grid.grid)
get_reffes(grid::GridWithFEMap) = get_reffes(grid.grid)
get_cell_type(grid::GridWithFEMap) = get_cell_type(grid.grid)

struct ModelWithFEMap{Dc,Dp,T} <: Gridap.Geometry.DiscreteModel{Dc,Dp}
  model::DiscreteModel{Dc,Dp}
  mapped_grid::GridWithFEMap{Dc,Dp}
  map::T
  function ModelWithFEMap(model::DiscreteModel{Dc,Dp},geo_map) where {Dc,Dp}
    model_map=get_cell_map(model)
    map=lazy_map(∘,geo_map,model_map)
    mapped_grid = GridWithFEMap(get_grid(model),map)
    new{Dc,Dp,typeof(map)}(model,mapped_grid,map)
  end
end

function ModelWithFEMap(model::DiscreteModel{Dc,Dp},uh::FEFunction) where {Dc,Dp}
  ModelWithFEMap(model,Gridap.CellData.get_data(uh))
end

import Gridap.Geometry.get_cell_map
import Gridap.Geometry.get_grid
import Gridap.Geometry.get_grid_topology
import Gridap.Geometry.get_face_labeling

get_cell_map(model::ModelWithFEMap) = model.geo_map
get_grid(model::ModelWithFEMap) = model.mapped_grid
get_grid_topology(model::ModelWithFEMap) = Gridap.Geometry.get_grid_topology(model.model)
get_face_labeling(model::ModelWithFEMap) = Gridap.Geometry.get_face_labeling(model.model)

function Grid(::Type{ReferenceFE{d}},model::ModelWithFEMap) where {d}
  get_grid(model)
end

# Basic tests

domain = (0,1)
partition = (2)
model = CartesianDiscreteModel(domain,partition)

T = VectorValue{1,Float64}; order = 1; pol = Polytope(HEX_AXIS...)
reffe = LagrangianRefFE(T,pol,order)

Vₕ = FESpace(model,reffe;conformity=:H1)
u(x) = 2*x            # Analytical solution (for Dirichlet data)
Uₕ = TrialFESpace(Vₕ,u)
uh = interpolate_everywhere(u,Uₕ)

geo_map = uh

new_model = ModelWithFEMap(model,uh)

grid = get_grid(model)
GridWithFEMap(grid,uh)

grid = get_grid(new_model)

get_node_coordinates(grid.grid)
get_node_coordinates(grid)

Tₕ = Triangulation(new_model)
get_node_coordinates(Tₕ)

test_discrete_model(new_model)
test_grid(grid)
test_triangulation(Tₕ)

# FE problem test

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

domain = (0,2)
partition = (2)
model = CartesianDiscreteModel(domain,partition)
uh = fe_problem(model)
get_free_dof_values(uh)

function get_map_model()
  domain = (0,1)
  partition = (2)
  model = CartesianDiscreteModel(domain,partition)

  T = VectorValue{1,Float64}; order = 1; pol = Polytope(HEX_AXIS...)
  reffe = LagrangianRefFE(T,pol,order)

  Vₕ = FESpace(model,reffe;conformity=:H1)
  d(x) = 2*x            # Analytical solution (for Dirichlet data)
  Uₕ = TrialFESpace(Vₕ,d)
  dh = interpolate_everywhere(d,Uₕ)
  new_model = ModelWithFEMap(model,dh)
end

map_model = get_map_model()
uh_map = fe_problem(map_model)

get_free_dof_values(uh_map)

end # module
