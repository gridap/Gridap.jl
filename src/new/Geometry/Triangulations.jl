
"""
"""
abstract type Triangulation{Dc,Dp} end

"""
"""
function get_cell_coordinates(trian::Triangulation)
  @abstractmethod
end

"""
"""
function get_reffes(trian::Triangulation)
  @abstractmethod
end

"""
"""
function get_cell_types(trian::Triangulation)
  @abstractmethod
end

"""
"""
function test_triangulation(trian::Triangulation{Dc,Dp}) where {Dc,Dp}
  @test num_cell_dims(trian) == Dc
  @test num_point_dims(trian) == Dp
  @test num_cell_dims(typeof(trian)) == Dc
  @test num_point_dims(typeof(trian)) == Dp
  cell_coords = get_cell_coordinates(trian)
  @test isa(cell_coords,AbstractArray{<:AbstractVector{<:Point}})
  reffes = get_reffes(trian)
  @test isa(reffes,AbstractVector{<:LagrangianRefFE})
  cell_types = get_cell_types(trian)
  @test isa(cell_types,AbstractArray{<:Integer})
  ncells = num_cells(trian)
  @test ncells == length(cell_coords)
  @test ncells == length(cell_types)
end

# Some API

"""
"""
num_cells(trian::Triangulation) = length(get_cell_types(trian))

"""
"""
num_cell_dims(::Triangulation{Dc,Dp}) where {Dc,Dp} = Dc
num_cell_dims(::Type{<:Triangulation{Dc,Dp}}) where {Dc,Dp} = Dc

"""
"""
num_point_dims(::Triangulation{Dc,Dp}) where {Dc,Dp} = Dp
num_point_dims(::Type{<:Triangulation{Dc,Dp}}) where {Dc,Dp} = Dp

"""

It is not desirable to iterate over the resulting array
for large number of cells. (The `LagrangianRefFE` type is
unstable deliberately to simplify its type signature)
"""
function get_cell_reffes(trian::Triangulation)
  type_to_reffe = get_reffes(trian)
  cell_to_type = get_cell_types(trian)
  _get_cell_data(type_to_reffe,cell_to_type)
end

function get_cell_shapefuns(trian::Triangulation)
  type_to_reffes = get_reffes(trian)
  cell_to_type = get_cell_types(trian)
  type_to_shapefuns = map(get_shapefuns, type_to_reffes)
  _get_cell_data(type_to_shapefuns,cell_to_type)
end

"""
"""
function get_cell_map(trian::Triangulation)
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_shapefuns = get_cell_shapefuns(trian)
  lincomb(cell_to_shapefuns, cell_to_coords)
end

# Helpers for Triangulation

function _get_cell_data(type_to_data, cell_to_type)
  CompressedArray(type_to_data,cell_to_type)
end

function _get_cell_data(type_to_data, cell_to_type::Fill)
  ncells = length(cell_to_type)
  @assert length(type_to_data) == 1 "Only one reference element expected"
  @assert cell_to_type.value == 1 "Only one type of reference element expected"
  data = first(type_to_data)
  Fill(data,ncells)
end

# ConformingTrian interface

"""
"""
abstract type ConformingTrian{Dc,Dp} <: Triangulation{Dc,Dp} end

"""
"""
function get_node_coordinates(trian::ConformingTrian)
  @abstractmethod
end

"""
"""
function get_cell_nodes(trian::ConformingTrian)
  @abstractmethod
end

# Methods from triangulation

function get_cell_coordinates(trian::ConformingTrian)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_nodes(trian)
  LocalToGlobalArray(cell_to_nodes,node_to_coords)
end

# Some API

num_nodes(trian::ConformingTrian) = length(get_node_coordinates(trian))

