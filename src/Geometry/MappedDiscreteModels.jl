function _cell_vector_to_dof_vector!(dof_vector,cell_node_ids, cell_vector)
  cache_cell_node_ids = array_cache(cell_node_ids)
  cache_cell_vector   = array_cache(cell_vector)
  for k=1:length(cell_node_ids)
    current_node_ids = getindex!(cache_cell_node_ids,cell_node_ids,k)
    current_values   = getindex!(cache_cell_vector,cell_vector,k)
    for (i,id) in enumerate(current_node_ids)
      dof_vector[current_node_ids[i]]=current_values[i]
    end
  end
end

# """
#     MappedGrid

# Represent a grid with a geometrical map that is the composition of
# a reference to a physical space map (standard)
# and a (vector-valued) map from physical space to physical space. E.g.,
# it can be used to include a high order map implemented by any map that is
# a `CellField`.
# """
"""
    struct MappedGrid{Dc,Dp,T,M,L} <: Grid{Dc,Dp}
"""
struct MappedGrid{Dc,Dp,T,M,L} <: Grid{Dc,Dp}
  grid::Grid{Dc,Dp}
  geo_map::T  # Composition of old map and new one
  phys_map::M # New map in the physical space
  node_coords::L
  function MappedGrid(grid::Grid{Dc,Dp},phys_map) where {Dc,Dp}

    @assert length(phys_map) == num_cells(grid)
    @assert eltype(phys_map) <: Field

    function _compute_node_coordinates(grid,phys_map)
      cell_node_ids = get_cell_node_ids(grid)
      old_nodes = get_node_coordinates(grid)
      node_coordinates = Vector{eltype(old_nodes)}(undef,length(old_nodes))
      c_coor = get_cell_coordinates(grid)
      map_c_coor = lazy_map(evaluate,phys_map,c_coor)
      _cell_vector_to_dof_vector!(node_coordinates,cell_node_ids,map_c_coor)
      return node_coordinates
    end

    model_map=get_cell_map(grid)
    geo_map=lazy_map(∘,phys_map,model_map)
    node_coords = collect(_compute_node_coordinates(grid,phys_map))
    new{Dc,Dp,typeof(geo_map),typeof(phys_map),typeof(node_coords)}(grid,geo_map,phys_map,node_coords)
  end
end

function MappedGrid(grid::Grid{Dc,Dp},phys_map::Function) where {Dc,Dp}
  c_map = Fill(GenericField(phys_map),num_cells(grid))
  MappedGrid(grid,c_map)
end

get_node_coordinates(grid::MappedGrid) = grid.node_coords
get_cell_node_ids(grid::MappedGrid) = get_cell_node_ids(grid.grid)
get_reffes(grid::MappedGrid) = get_reffes(grid.grid)
get_cell_type(grid::MappedGrid) = get_cell_type(grid.grid)

"""
  struct MappedDiscreteModel{Dc,Dp} <: DiscreteModel{Dc,Dp}

Represent a model with a `MappedGrid` grid.
See also [`MappedGrid`](@ref).
"""
struct MappedDiscreteModel{Dc,Dp} <: DiscreteModel{Dc,Dp}
  model::DiscreteModel{Dc,Dp}
  mapped_grid::MappedGrid{Dc,Dp}
  function MappedDiscreteModel(model::DiscreteModel{Dc,Dp},phys_map) where {Dc,Dp}
    mapped_grid = MappedGrid(get_grid(model),phys_map)
    new{Dc,Dp}(model,mapped_grid)
  end
end

get_grid(model::MappedDiscreteModel) = model.mapped_grid
get_cell_map(model::MappedDiscreteModel) = get_cell_map(model.mapped_grid)
get_grid_topology(model::MappedDiscreteModel) = get_grid_topology(model.model)
get_face_labeling(model::MappedDiscreteModel) = get_face_labeling(model.model)

function Grid(::Type{ReferenceFE{d}},model::MappedDiscreteModel) where {d}
  get_grid(model)
end
