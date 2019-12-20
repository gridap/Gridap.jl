
"""
    abstract type Triangulation{Dc,Dp}

Abstract type representing an arbitrary tiling, tessellation,
or triangulation of a domain of parametric dimension `Dc` and
physical dimension `Dp`.

We define a triangulation from two basic ingredients: 

- the cell-wise nodal coordinates of the cells in the triangulation, plus
- an interpolation of this cell-wise coordinates into the cells interior.

Note that this type represents general triangulations (not necessarily conforming),
which is the minimum geometrical information needed to perform cell-wise numerical integration.

The `Triangulation` interface is defined by overloading these methods:


- [`get_cell_coordinates(trian::Triangulation)`](@ref)
- [`get_reffes(trian::Triangulation)`](@ref)
- [`get_cell_type(trian::Triangulation)`](@ref)

and it can be tested with

- [`test_triangulation`](@ref)

"""
abstract type Triangulation{Dc,Dp} end

"""
    get_cell_coordinates(trian::Triangulation) -> AbstractArray{Vector{<:Point{Dp}}}
"""
function get_cell_coordinates(trian::Triangulation)
  @abstractmethod
end

"""
    get_reffes(trian::Triangulation) -> Vector{LagrangianRefFE}
"""
function get_reffes(trian::Triangulation)
  @abstractmethod
end

"""
    get_cell_type(trian::Triangulation) -> AbstractVector{<:Integer}
"""
function get_cell_type(trian::Triangulation)
  @abstractmethod
end

"""
    test_triangulation(trian::Triangulation)
"""
function test_triangulation(trian::Triangulation{Dc,Dp}) where {Dc,Dp}
  @test num_cell_dims(trian) == Dc
  @test num_point_dims(trian) == Dp
  @test num_cell_dims(typeof(trian)) == Dc
  @test num_point_dims(typeof(trian)) == Dp
  cell_coords = get_cell_coordinates(trian)
  @test isa(cell_coords,AbstractArray{<:AbstractVector{<:Point}})
  reffes = get_reffes(trian)
  @test isa(reffes,AbstractVector{LagrangianRefFE{Dc}})
  cell_types = get_cell_type(trian)
  @test isa(cell_types,AbstractArray{<:Integer})
  ncells = num_cells(trian)
  @test ncells == length(cell_coords)
  @test ncells == length(cell_types)
end

# Some API

"""
    num_cells(trian::Triangulation) -> Int
"""
num_cells(trian::Triangulation) = length(get_cell_type(trian))

"""
    num_cell_dims(::Triangulation) -> Int
    num_cell_dims(::Type{<:Triangulation}) -> Int
"""
num_cell_dims(::Triangulation{Dc,Dp}) where {Dc,Dp} = Dc
num_cell_dims(::Type{<:Triangulation{Dc,Dp}}) where {Dc,Dp} = Dc

"""
    num_point_dims(::Triangulation) -> Int
    num_point_dims(::Type{<:Triangulation}) -> Int
"""
num_point_dims(::Triangulation{Dc,Dp}) where {Dc,Dp} = Dp
num_point_dims(::Type{<:Triangulation{Dc,Dp}}) where {Dc,Dp} = Dp

"""
    num_dims(::Triangulation) -> Int
    num_dims(::Type{<:Triangulation}) -> Int

Equivalent to `num_cell_dims`.
"""
num_dims(g::Triangulation{Dc}) where Dc = Dc
num_dims(::Type{<:Triangulation{Dc}}) where Dc = Dc

"""
    is_affine(trian::Triangulation) -> Bool
"""
function is_affine(trian::Triangulation)
  reffes = get_reffes(trian)
  all(map(is_affine,reffes))
end

"""
    is_first_order(trian::Triangulation) -> Bool
"""
function is_first_order(trian::Triangulation)
  reffes = get_reffes(trian)
  all(map(is_first_order,reffes))
end

"""
    get_cell_reffes(trian::Triangulation) -> Vector{<:NodalReferenceFEs}

It is not desirable to iterate over the resulting array
for large number of cells if the underlying reference FEs
are of different Julia type.
"""
function get_cell_reffes(trian::Triangulation)
  type_to_reffe = get_reffes(trian)
  cell_to_type = get_cell_type(trian)
  _get_cell_data(type_to_reffe,cell_to_type)
end

"""
    get_cell_shapefuns(trian::Triangulation) -> Vector{<:Field}
"""
function get_cell_shapefuns(trian::Triangulation)
  type_to_reffes = get_reffes(trian)
  cell_to_type = get_cell_type(trian)
  type_to_shapefuns = map(get_shapefuns, type_to_reffes)
  _get_cell_data(type_to_shapefuns,cell_to_type)
end

"""
    get_cell_map(trian::Triangulation) -> Vector{<:Field}
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

