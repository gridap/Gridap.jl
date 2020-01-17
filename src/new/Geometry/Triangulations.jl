
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

Optional interface:

For triangulations living in a space of co-dimension 1, the following method can be defined:

- [`get_normal_vector(trian::Triangulation)`]

In some cases, concrete implementations want to override the default implementation of the following methods:

- [`restrict(f::AbstractArray, trian::Triangulation)`]
- [`get_cell_id(f::AbstractArray, trian::Triangulation)`]

The (mandatory) `Triangulation` interface can be tested with

- [`test_triangulation`](@ref)


"""
abstract type Triangulation{Dc,Dp} <: GridapType end

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
    get_normal_vector(trian::Triangulation)
"""
function get_normal_vector(trian::Triangulation{Dc,Dp}) where {Dc,Dp}
  if Dp == Dc + 1
    @abstractmethod
  else
    @unreachable "get_normal_vector does not make sense for this triangulation"
  end
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
    reindex(a::AbstractArray, trian::Triangulation)
"""
function reindex(a::AbstractArray,trian::Triangulation)
  reindex(a,get_cell_id(trian))
end

"""
    get_cell_id(trian::Triangulation)
"""
function get_cell_id(trian::Triangulation)
  identity_vector(num_cells(trian))
end

"""
    restrict(f::AbstractArray, trian::Triangulation)
"""
function restrict(f::AbstractArray,trian::Triangulation)
  f
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

"""
"""
function CellField(object,trian::Triangulation)
  cm = get_cell_map(trian)
  convert_to_cell_field(object,cm)
end

"""
In contrast to get_cell_map, the returned object:
- is a CellField
- its gradient is the identity tensor
"""
function get_physical_coordinate(trian::Triangulation)
  CellField(_phys_coord,trian)
end

_phys_coord(x) = x

_phys_coord_grad(x) = one(typeof(outer(x,x)))

gradient(::typeof(_phys_coord)) = _phys_coord_grad

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

"""
"""
function restrict(cf::CellField,trian::Triangulation)
  a = get_array(cf)
  r = restrict(a,trian)
  _restrict_cell_field(r,trian)
end

"""
"""
struct SkeletonCellField{L,R}
  left::L
  right::R
end

function jump(sf::SkeletonCellField)
  sf.left - sf.right
end

function mean(sf::SkeletonCellField)
  operate_cell_field(_mean,sf.left,sf.right)
end

_mean(x,y) = 0.5*x + 0.5*y

function _restrict_cell_field(r::SkeletonPair,trian)
  cm = get_cell_map(trian)
  la = r.left
  ra = r.right
  l = GenericCellField(la,cm)
  r = GenericCellField(ra,cm)
  SkeletonCellField(l,r)
end

function _restrict_cell_field(r::AbstractArray,trian)
  cm = get_cell_map(trian)
  GenericCellField(r,cm)
end

