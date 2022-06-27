"""
Given a Discrete Model and a reffe, builds a new grid in which the geometrical
map is a `FEFunction`. This is useful when considering geometrical maps that are
the result of a FE problem (mesh displacement).
"""
struct GridWithFEMap{Dc,Dp,A,B,C,D,E,F} <: Grid{Dc,Dp}
  grid::Grid{Dc,Dp}
  fe_sp::A
  scal_fe_sp::B
  fe_map::C
  node_coords::D
  reffes::E
  cell_type::F
end

function GridWithFEMap(model,reffe::ReferenceFE; kwargs...)
  cell_type = Fill(1,num_cells(model))
  _reffes = [reffe]
  reffes = lazy_map(Reindex(_reffes),cell_type)
  GridWithFEMap(model,cell_type,reffes; kwargs...)
end

function GridWithFEMap(model,cell_type,reffes::AbstractArray{<:ReferenceFE}; kwargs...)

  # Create a FESpace for the geometrical description
  # The FEFunction that describes the coordinate field
  # or displacement can be a interpolation or a solution
  # of a mesh displacement problem
  Vₕ = FESpace(model,reffes;conformity=:H1,kwargs...)

  os = lazy_map(get_order,reffes)
  ps = lazy_map(get_polytope,reffes)
  f(a,b) = LagrangianRefFE(Float64,a,b)
  s_reffes = lazy_map(f,ps,os)
  Vₕ_scal = FESpace(model,s_reffes;conformity=:H1)

  grid = get_grid(model)
  geo_map = get_cell_map(grid)

  cell_ctype = get_cell_type(grid)
  c_reffes = get_reffes(grid)

  # Create a fe_map using the cell_map that can be evaluated at the
  # vertices of the fe space (nodal type)
  # This returns a FEFunction initialised with the coordinates
  # But this is to build the FEFunction that will be inserted, it is
  # an advanced constructor, not needed at this stage

  c_dofs = get_fe_dof_basis(Vₕ)
  dof_basis = get_data(c_dofs)
  c_nodes = lazy_map(get_nodes,get_data(c_dofs))

  xh = zero(Vₕ)
  c_dofv = lazy_map(evaluate,dof_basis,geo_map)

  Uₕ = TrialFESpace(Vₕ)

  fv = get_free_dof_values(xh)
  dv = get_dirichlet_dof_values(get_fe_space(xh))
  gather_free_and_dirichlet_values!(fv,dv,Uₕ,c_dofv)

  c_xh = lazy_map(evaluate,get_data(xh),c_nodes)
  c_scal_ids = get_cell_dof_ids(Vₕ_scal)


  nodes_coords = Vector{eltype(eltype(c_xh))}(undef,num_free_dofs(Vₕ_scal))
  Geometry._cell_vector_to_dof_vector!(nodes_coords,c_scal_ids,c_xh)

  GridWithFEMap(grid, Vₕ, Vₕ_scal, xh, nodes_coords, reffes, cell_type)
end

function _compute_node_coordinates(grid,xh)
  c_dofs = get_fe_dof_basis(grid.fe_sp)
  c_nodes = lazy_map(get_nodes,get_data(c_dofs))

  c_xh = lazy_map(evaluate,get_data(xh),c_nodes)
  c_scal_ids = get_cell_dof_ids(grid.scal_fe_sp)

  nodes_coords = grid.node_coords
  Geometry._cell_vector_to_dof_vector!(nodes_coords,c_scal_ids,c_xh)

  return nodes_coords
end

function add_mesh_displacement!(grid::GridWithFEMap,dh::FEFunction)

  Xh = grid.fe_map

  Xh_fv = get_free_dof_values(Xh)
  Xh_dv = get_dirichlet_dof_values(get_fe_space(Xh))

  Xh_fv .= Xh_fv .+ get_free_dof_values(dh)
  Xh_dv .= Xh_dv .+ get_dirichlet_dof_values(get_fe_space(dh))

  _compute_node_coordinates(grid,Xh)
end

function update_coordinates!(grid::GridWithFEMap,dh::FEFunction)

  Xh = grid.fe_map

  Xh_fv = get_free_dof_values(Xh)
  Xh_dv = get_dirichlet_dof_values(get_fe_space(Xh))

  Xh_fv .= get_free_dof_values(dh)
  Xh_dv .= get_dirichlet_dof_values(get_fe_space(dh))

  _compute_node_coordinates(grid,Xh)

end

# Interface
get_node_coordinates(grid::GridWithFEMap) = grid.node_coords
get_cell_node_ids(grid::GridWithFEMap) = get_cell_dof_ids(grid.scal_fe_sp)
get_reffes(grid::GridWithFEMap) = grid.reffes
get_cell_type(grid::GridWithFEMap) = grid.cell_type
OrientationStyle(grid::GridWithFEMap) = OrientationStyle(grid.grid)
RegularityStyle(grid::GridWithFEMap) = RegularityStyle(grid.grid)
# santiagobadia: To think (does it make sense here, I think it is for surface problems)
get_facet_normal(grid::GridWithFEMap) = get_facet_normal(grid.grid)


struct DiscreteModelWithFEMap{Dc,Dp} <: DiscreteModel{Dc,Dp}
  model::DiscreteModel{Dc,Dp}
  mapped_grid::GridWithFEMap{Dc,Dp}
  function DiscreteModelWithFEMap(model::DiscreteModel{Dc,Dp},phys_map) where {Dc,Dp}
    mapped_grid = GridWithFEMap(get_grid(model),phys_map)
    new{Dc,Dp}(model,mapped_grid)
  end
end

function DiscreteModelWithFEMap(model::DiscreteModel,reffe::ReferenceFE; kwargs...)
  mapped_grid = GridWithFEMap(model,reffe;kwargs...)
  MappedDiscreteModel(model,mapped_grid)
end

function DiscreteModelWithFEMap(model,cell_type,reffes::AbstractArray{<:ReferenceFE}; kwargs...)
  mapped_grid = GridWithFEMap(model,cell_type,reffes;kwargs...)
  MappedDiscreteModel(model,mapped_grid)
end
