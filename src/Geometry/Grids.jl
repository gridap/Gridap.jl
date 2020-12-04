
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

abstract type OrientationStyle end
struct Oriented <: OrientationStyle end
struct NonOriented <: OrientationStyle end

"""
    OrientationStyle(::Type{<:Grid})
    OrientationStyle(::Grid)

`Oriented()` if has oriented faces, `NonOriented()` otherwise (default).
"""
OrientationStyle(a::Grid) = OrientationStyle(typeof(a))
OrientationStyle(::Type{<:Grid}) = NonOriented()

abstract type RegularityStyle end
struct Regular <: RegularityStyle end
struct Irregular <: RegularityStyle end

"""
    RegularityStyle(::Type{<:Grid})
    RegularityStyle(::Grid)

`Regular()` if no hanging-nodes default), `Irregular()` otherwise.
"""
RegularityStyle(a::Grid) = RegularityStyle(typeof(a))
RegularityStyle(::Type{<:Grid}) = Regular()

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
  @test OrientationStyle(trian) in (Oriented(), NonOriented())
  @test RegularityStyle(trian) in (Regular(), Irregular())
  @test is_oriented(trian) == (OrientationStyle(trian) == Oriented())
  @test is_regular(trian) == (RegularityStyle(trian) == Regular())
end

# Methods from triangulation

function get_cell_coordinates(trian::Grid)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_nodes(trian)
  lazy_map(Broadcasting(Reindex(node_to_coords)),cell_to_nodes)
end

# Some API

"""
    is_oriented(::Type{<:Grid}) -> Bool
    is_oriented(a::Grid) -> Bool
"""
is_oriented(a::Grid) = is_oriented(typeof(a))
is_oriented(a::Type{T}) where T<:Grid = OrientationStyle(T) == Oriented()

"""
    is_regular(::Type{<:Grid}) -> Bool
    is_regular(a::Grid) -> Bool
"""
is_regular(a::Grid) = is_regular(typeof(a))
is_regular(a::Type{T}) where T<:Grid = RegularityStyle(T) == Regular()

"""
    Grid(reffe::LagrangianRefFE)
"""
function Grid(reffe::LagrangianRefFE)
  UnstructuredGrid(reffe)
end

"""
    compute_linear_grid(reffe::LagrangianRefFE)
"""
function compute_linear_grid(reffe::LagrangianRefFE)
  @notimplementedif ! is_n_cube(get_polytope(reffe)) "linear grid only implemented for n-cubes at the moment"
  if get_order(reffe) == 0
    D = num_cell_dims(reffe)
    partition = tfill(1,Val{D}())
  else
    partition = get_orders(reffe)
  end
  pmin, pmax = get_bounding_box(get_polytope(reffe))
  desc = CartesianDescriptor(pmin,pmax,partition)
  CartesianGrid(desc)
end

"""
    compute_reference_grid(p::LagrangianRefFE, nelems::Integer)
"""
function compute_reference_grid(reffe::LagrangianRefFE, nelems::Integer)
  p = get_polytope(reffe)
  r = LagrangianRefFE(Float64,p,nelems)
  compute_linear_grid(r)
end

"""
    Grid(::Type{ReferenceFE{d}},p::Polytope) where d
"""
function Grid(::Type{ReferenceFE{d}},p::Polytope) where d
  UnstructuredGrid(ReferenceFE{d},p)
end

"""
    simplexify(grid::Grid)
"""
function simplexify(grid::Grid)
  simplexify(UnstructuredGrid(grid))
end

