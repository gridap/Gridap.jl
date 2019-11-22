
"""
"""
abstract type ConformingTriangulation{Dc,Dp} <: Triangulation{Dc,Dp} end

"""
"""
function get_node_coordinates(trian::ConformingTriangulation)
  @abstractmethod
end

"""
"""
function get_cell_nodes(trian::ConformingTriangulation)
  @abstractmethod
end

"""
"""
function test_conforming_triangulation(trian::ConformingTriangulation)
  test_triangulation(trian)
  nodes_coords = get_node_coordinates(trian)
  @test isa(nodes_coords,AbstractArray{<:Point})
  cell_nodes = get_cell_nodes(trian)
  @test isa(cell_nodes,AbstractArray{<:AbstractArray{<:Integer}})
  @test num_nodes(trian) == length(nodes_coords)
end

# Methods from triangulation

function get_cell_coordinates(trian::ConformingTriangulation)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_nodes(trian)
  LocalToGlobalArray(cell_to_nodes,node_to_coords)
end

# Some API

num_nodes(trian::ConformingTriangulation) = length(get_node_coordinates(trian))

