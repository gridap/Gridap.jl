
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

# Mock triangulation for testing purposes

struct ConformingTrianMock <: ConformingTriangulation{2,2} end

function get_node_coordinates(::ConformingTrianMock)
  Point{2,Float64}[(0,0),(1,0),(2,0),(1,1),(2,1),(0,2),(2,2)]
end

function get_cell_nodes(::ConformingTrianMock)
  [[1,2,6,4],[2,3,4,5],[4,5,7],[4,7,6]]
end

function get_reffes(::ConformingTrianMock)
  order = 1
  tri3 = LagrangianRefFE(Float64,TRI,order)
  quad4 = LagrangianRefFE(Float64,QUAD,order)
  [tri3, quad4]
end

function get_cell_types(::ConformingTrianMock)
  [2,2,1,1]
end

