
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
- [`ConformityStyle(::Type{<:Grid})`](@ref)

The interface of `Grid` is tested with

- [`test_conforming_triangulation`](@ref)

"""
abstract type Grid{Dc,Dp} <: Triangulation{Dc,Dp} end


"""
    is_oriented(::Type{<:Grid}) -> Bool
    is_oriented(a::Grid) -> Bool
"""
is_oriented(a::Grid) = _is_oriented(OrientationStyle(a))
is_oriented(a::Type{<:Grid}) = _is_oriented(OrientationStyle(a))
_is_oriented(::Val{true}) = true
_is_oriented(::Val{false}) = false

"""
    abstract type ConformityStyle end
"""
abstract type ConformityStyle end

"""
    struct RegularConformity <: ConformityStyle end
"""
struct RegularConformity <: ConformityStyle end

"""
    struct IrregularHConformity <: ConformityStyle end
"""
struct IrregularHConformity <: ConformityStyle end

"""
    struct IrregularPConformity <: ConformityStyle end
"""
struct IrregularPConformity <: ConformityStyle end

"""
    struct IrregularHPConformity <: ConformityStyle end
"""
struct IrregularHPConformity <: ConformityStyle end

"""
    ConformityStyle(::Type{<:Grid})
    ConformityStyle(a::Grid)
"""
ConformityStyle(::Type{<:Grid}) = RegularConformity()
ConformityStyle(a::Grid) = ConformityStyle(typeof(a))

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
    test_conforming_triangulation(trian::Grid)
"""
function test_conforming_triangulation(trian::Grid)
  test_triangulation(trian)
  nodes_coords = get_node_coordinates(trian)
  @test isa(nodes_coords,AbstractArray{<:Point})
  cell_nodes = get_cell_nodes(trian)
  @test isa(cell_nodes,AbstractArray{<:AbstractArray{<:Integer}})
  @test num_nodes(trian) == length(nodes_coords)
  @test isa(is_oriented(trian),Bool)
  @test OrientationStyle(trian) in (Val{false}(), Val{true}())
  @test isa(ConformityStyle(trian),ConformityStyle)
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
    Grid(reffe::NodalReferenceFE)
"""
function Grid(reffe::NodalReferenceFE)
  UnstructuredGrid(reffe)
end

"""
    Grid(::Type{<:ReferenceFE{d}},p::Polytope) where d
"""
function Grid(::Type{<:ReferenceFE{d}},p::Polytope) where d
  UnstructuredGrid(NodalReferenceFE{d},p)
end

"""
    Grid(::Type{<:ReferenceFE{d}},trian::Grid) where d
"""
function Grid(::Type{<:ReferenceFE{d}},trian::Grid) where d
  UnstructuredGrid(NodalReferenceFE{d},trian)
end

"""
    replace_reffes(grid::Grid,reffes::Vector{<:NodalReferenceFE})
"""
function replace_reffes(grid::Grid,reffes::Vector{<:NodalReferenceFE})
  model = UnstructuredDiscreteModel(grid)
  model2 = replace_reffes(model,reffes)
  D = num_cell_dims(grid)
  Grid(ReferenceFE{D},model2)
end

