# module ModelWithFEMaps

# using Gridap
# using Gridap.Arrays
# using Gridap.Fields
# using Gridap.ReferenceFEs
# using Gridap.Geometry
# using Gridap.CellData
# using Gridap.FESpaces

# import Gridap.ReferenceFEs.get_node_coordinates
# import Gridap.Geometry.get_cell_node_ids
# import Gridap.Geometry.get_reffes
# import Gridap.Geometry.get_cell_type
# import Gridap.Geometry.Grid

# import Gridap.Geometry.get_cell_map
# import Gridap.Geometry.get_grid
# import Gridap.Geometry.get_grid_topology
# import Gridap.Geometry.get_face_labeling

# """
#     MappedGrid

# Reprepresent a grid with a geometrical map that is the composition of
# a reference to a physical space map (standard)
# and a (vector-valued) map from physical space to physical space. E.g.,
# it can be used to include a high order map implemented by any map that is
# a `CellField`.
# """
struct MappedGrid{Dc,Dp,T} <: Grid{Dc,Dp}
  grid::Grid{Dc,Dp}
  geo_map::T
  function MappedGrid(grid::Grid{Dc,Dp},geo_map) where {Dc,Dp}
    model_map=get_cell_map(grid)
    map=lazy_map(∘,geo_map,model_map)
    new{Dc,Dp,typeof(map)}(grid,map)
  end
end

function get_node_coordinates(grid::MappedGrid)
  parent_node_to_coords = get_node_coordinates(grid.grid)
  lazy_map(grid.geo_map,parent_node_to_coords)
end

get_cell_node_ids(grid::MappedGrid) = get_cell_node_ids(grid.grid)
get_reffes(grid::MappedGrid) = get_reffes(grid.grid)
get_cell_type(grid::MappedGrid) = get_cell_type(grid.grid)

"""
    MappedDiscreteModel

Reprepresent a model with a `MappedGrid` grid.
See also [`MappedGrid`](@ref).
"""
struct MappedDiscreteModel{Dc,Dp,T} <: DiscreteModel{Dc,Dp}
  model::DiscreteModel{Dc,Dp}
  mapped_grid::MappedGrid{Dc,Dp}
  map::T
  function MappedDiscreteModel(model::DiscreteModel{Dc,Dp},geo_map) where {Dc,Dp}
    model_map=get_cell_map(model)
    map=lazy_map(∘,geo_map,model_map)
    mapped_grid = MappedGrid(get_grid(model),map)
    new{Dc,Dp,typeof(map)}(model,mapped_grid,map)
  end
end


get_cell_map(model::MappedDiscreteModel) = model.geo_map
get_grid(model::MappedDiscreteModel) = model.mapped_grid
get_grid_topology(model::MappedDiscreteModel) = get_grid_topology(model.model)
get_face_labeling(model::MappedDiscreteModel) = get_face_labeling(model.model)

function Grid(::Type{ReferenceFE{d}},model::MappedDiscreteModel) where {d}
  get_grid(model)
end

# end # module
