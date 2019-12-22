
"""
    abstract type Grid{Dc,Dp} <: Triangulation{Dc,Dp}

Abstract type that represents conforming triangulations, whose cell-wise
nodal coordinates are defined with a vector of nodal coordinates, plus
a cell-wise vector of node ids.

The interface of `Grid` is defined by overloading the
methods in `Triangulation` plus the following ones:

- [`get_node_coordinates(trian::Grid)`](@ref)
- [`get_cell_nodes(trian::Grid)`](@ref)

From these two methods a default implementation of [`get_cell_coordinates(trian::Triangulation)`](@ref)
is available.

The `Grid`  interface has the following traits

- [`OrientationStyle(::Type{<:Grid})`](@ref)
- [`RegularityStyle(::Type{<:Grid})`](@ref)

The interface of `Grid` is tested with

- [`test_grid`](@ref)

"""
abstract type Grid{Dc,Dp} <: Triangulation{Dc,Dp} end

# Traits

"""
    OrientationStyle(::Type{<:Grid}) -> Val{Bool}
    OrientationStyle(::Grid) -> Val{Bool}

`Val{true}()` if has oriented faces, `Val{false}()` otherwise (default).
"""
OrientationStyle(::Type{<:Grid}) = Val{false}()
OrientationStyle(a::Grid) = OrientationStyle(typeof(a))

"""
    RegularityStyle(::Type{<:Grid}) -> Val{Bool}
    RegularityStyle(::Grid) -> Val{Bool}

`Val{true}()` if no hanging-nodes (refault), `Val{false}()` otherwise.
"""
RegularityStyle(::Type{<:Grid}) = Val{true}()
RegularityStyle(a::Grid) = RegularityStyle(typeof(a))

# Interface

"""
    get_node_coordinates(trian::Grid) -> AbstractArray{<:Point{Dp}}
"""
function get_node_coordinates(trian::Grid)
  @abstractmethod
end

"""
    get_cell_nodes(trian::Grid)
"""
function get_cell_nodes(trian::Grid)
  @abstractmethod
end

"""
    test_grid(trian::Grid)
"""
function test_grid(trian::Grid)
  test_triangulation(trian)
  nodes_coords = get_node_coordinates(trian)
  @test isa(nodes_coords,AbstractArray{<:Point})
  cell_nodes = get_cell_nodes(trian)
  @test isa(cell_nodes,AbstractArray{<:AbstractArray{<:Integer}})
  @test num_nodes(trian) == length(nodes_coords)
  @test isa(is_oriented(trian),Bool)
  @test isa(is_regular(trian),Bool)
  @test OrientationStyle(trian) in (Val{false}(), Val{true}())
  @test RegularityStyle(trian) in (Val{false}(), Val{true}())
end

# Methods from triangulation

function get_cell_coordinates(trian::Grid)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_nodes(trian)
  LocalToGlobalArray(cell_to_nodes,node_to_coords)
end

# Some API

"""
    num_nodes(trian::Grid) -> Int
"""
num_nodes(trian::Grid) = length(get_node_coordinates(trian))

"""
    is_oriented(::Type{<:Grid}) -> Bool
    is_oriented(a::Grid) -> Bool
"""
is_oriented(a::Grid) = _is_oriented(OrientationStyle(a))
is_oriented(a::Type{<:Grid}) = _is_oriented(OrientationStyle(a))

"""
    is_regular(::Type{<:Grid}) -> Bool
    is_regular(a::Grid) -> Bool
"""
is_regular(a::Grid) = _is_regular(RegularityStyle(a))
is_regular(a::Type{<:Grid}) = _is_regular(RegularityStyle(a))

"""
    Grid(reffe::LagrangianRefFE)
"""
function Grid(reffe::LagrangianRefFE)
  UnstructuredGrid(reffe)
end

"""
    Grid(::Type{ReferenceFE{d}},p::Polytope) where d
"""
function Grid(::Type{ReferenceFE{d}},p::Polytope) where d
  UnstructuredGrid(ReferenceFE{d},p)
end

"""
    GridTopology(grid::Grid)
    GridTopology(grid::Grid, cell_to_vertices::Table, vertex_to_node::Vector)
"""
function GridTopology(grid::Grid)
  _grid = UnstructuredGrid(grid)
  UnstructuredGridTopology(_grid)
end

function GridTopology(grid::Grid, cell_to_vertices::Table, vertex_to_node::AbstractVector)
  _grid = UnstructuredGrid(grid)
  UnstructuredGridTopology(_grid,cell_to_vertices,vertex_to_node)
end


