
"""
    abstract type ConformingTriangulation{Dc,Dp} <: Triangulation{Dc,Dp}

Abstract type that represents conforming triangulations, whose cell-wise
nodal coordinates are defined with a vector of nodal coordinates, plus
a cell-wise vector of node ids.

The interface of `ConformingTriangulation` is defined by overloading the
methods in `Triangulation` plus the following ones:

- [`get_node_coordinates(trian::ConformingTriangulation)`](@ref)
- [`get_cell_nodes(trian::ConformingTriangulation)`](@ref)

From these two methods a default implementation of [`get_cell_coordinates(trian::Triangulation)`](@ref)
is available.

The interface of `ConformingTriangulation` is tested with

- [`test_conforming_triangulation`](@ref)

"""
abstract type ConformingTriangulation{Dc,Dp} <: Triangulation{Dc,Dp} end

"""
    get_node_coordinates(trian::ConformingTriangulation) -> AbstractArray{<:Point{Dp}}
"""
function get_node_coordinates(trian::ConformingTriangulation)
  @abstractmethod
end

"""
    get_cell_nodes(trian::ConformingTriangulation)
"""
function get_cell_nodes(trian::ConformingTriangulation)
  @abstractmethod
end

"""
    test_conforming_triangulation(trian::ConformingTriangulation)
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

"""
    num_nodes(trian::ConformingTriangulation) -> Int
"""
num_nodes(trian::ConformingTriangulation) = length(get_node_coordinates(trian))


"""
    ConformingTriangulation(reffe::NodalReferenceFE)
"""
function ConformingTriangulation(reffe::NodalReferenceFE)
  UnstructuredGrid(reffe)
end

"""
    ConformingTriangulation(::Type{<:ReferenceFE{d}},p::Polytope) where d
"""
function ConformingTriangulation(::Type{<:ReferenceFE{d}},p::Polytope) where d
  UnstructuredGrid(NodalReferenceFE{d},p)
end

"""
    ConformingTriangulation(::Type{<:ReferenceFE{d}},trian::ConformingTriangulation) where d
"""
function ConformingTriangulation(::Type{<:ReferenceFE{d}},trian::ConformingTriangulation) where d
  UnstructuredGrid(NodalReferenceFE{d},trian)
end

