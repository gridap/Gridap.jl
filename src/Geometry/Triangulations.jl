
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

- [`get_facet_normal(trian::Triangulation)`]

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
    get_facet_normal(trian::Triangulation)
"""
function get_facet_normal(trian::Triangulation)
  Dp = num_point_dims(trian)
  Dc = num_cell_dims(trian)
  if Dp == Dc + 1
    @abstractmethod
  else
    @unreachable "get_facet_normal does not make sense for this triangulation"
  end
end

"""
If the reference (resp. physical) space is the same for both triangulaitons.
"""
have_compatible_domains(a::Triangulation,b::Triangulation) = a===b

# Trait that signals if the triangulation is a sub-mesh of a background triangulation

abstract type TriangulationStyle end
struct BackgroundTriangulation <: TriangulationStyle end
struct SubTriangulation <: TriangulationStyle end

TriangulationStyle(::Type{<:Triangulation}) = BackgroundTriangulation()
TriangulationStyle(::T) where T<:Triangulation = TriangulationStyle(T)

# Default methods needed to be overloaded by Triangulations that are integration sub-meshes

"""
Returns the background triangulation.
"""
get_background_triangulation(trian::Triangulation) =
  get_background_triangulation(trian,TriangulationStyle(trian))
get_background_triangulation(trian::Triangulation,::BackgroundTriangulation) = trian
get_background_triangulation(trian::Triangulation,::SubTriangulation) = @abstractmethod

#"""
#    reindex(a::AbstractArray, trian::Triangulation)
#"""
#function reindex(a::AbstractArray,trian::Triangulation)
#  reindex(a,get_cell_id(trian))
#end

"""
    get_cell_id(trian::Triangulation)

Map from the indices in the sub-triangulation to the indices in the background triangulation
"""
get_cell_id(trian::Triangulation) = get_cell_id(trian,TriangulationStyle(trian))
get_cell_id(trian::Triangulation,::BackgroundTriangulation) = IdentityVector(num_cells(trian))
get_cell_id(trian::Triangulation,::SubTriangulation) = @abstractmethod

#"""
#    restrict(f::AbstractArray, trian::Triangulation)
#"""
#function restrict(f::AbstractArray,trian::Triangulation)
#  f
#end

"""
Return the cell-wise map that goes from the reference space of the sub-triangulation to
the reference space of the background triangulation
"""
get_cell_ref_map(trian::Triangulation) = get_cell_ref_map(trian,TriangulationStyle(trian))
get_cell_ref_map(trian::Triangulation,::BackgroundTriangulation) = Fill(GenericField(identity),num_cells(trian))
get_cell_ref_map(trian::Triangulation,::SubTriangulation) = @abstractmethod

#"""
#Given an array aligned with the cells in the background triangulation, return another array
#aligned with the cells of the sub-triangulation. Do nothing if `trian` is already a 
#background triangulation.
#"""
#change_cell_index(a::AbstractArray,trian::Triangulation) = change_cell_index(a,trian,TriangulationStyle(trian))
#function change_cell_index(a::AbstractArray,trian::Triangulation,::BackgroundTriangulation)
#  @assert length(a) == num_cells(trian)
#  a
#end
#function change_cell_index(a::AbstractArray,trian::Triangulation,::SubTriangulation)
#  bgtrian = get_background_triangulation(trian)
#  @assert length(a) == num_cells(bgtrian)
#  lazy_map(Reindex(a),get_cell_id(trian))
#end
#
#"""
#Given an array (of arrays) of Field objects that is aligned with the sub-triangulation but
#the domain space of those fields corresponds to the reference space of the background triangulation,
#return another array (of arrays) of Field objects with domain in the sub-triangulation.
#Do nothing if `trian` is already a background triangulation. 
#Note that one typically needs to call `change_cell_index` before calling `change_cell_domain`.
#"""
#change_cell_domain(a::AbstractArray,trian::Triangulation) = change_cell_domain(a,trian,TriangulationStyle(trian))
#function change_cell_domain(a::AbstractArray,trian::Triangulation,::BackgroundTriangulation)
#  @assert length(a) == num_cells(trian)
#  a
#end
#function change_cell_domain(a::AbstractArray,trian::Triangulation,::SubTriangulation)
#  @assert length(a) == num_cells(trian)
#  cell_ref_map = get_cell_ref_map(trian)
#  lazy_map(Broadcasting(∘),a,cell_map)
#end

"""
    struct SkeletonPair{L,R} <: GridapType
      plus::L
      minus::R
    end
"""
struct SkeletonPair{L,R} <: GridapType
  plus::L
  minus::R
end

function Base.getproperty(x::SkeletonPair, sym::Symbol)
  if sym == :⁺
    x.plus
  elseif sym == :⁻
    x.minus
  else
    getfield(x, sym)
  end
end

function Base.propertynames(x::SkeletonPair, private=false)
  (fieldnames(typeof(x))...,:⁺,:⁻)
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
  @test isa(reffes,AbstractVector{<:LagrangianRefFE{Dc}})
  cell_types = get_cell_type(trian)
  @test isa(cell_types,AbstractArray{<:Integer})
  ncells = num_cells(trian)
  @test ncells == length(cell_coords)
  @test ncells == length(cell_types)
  @test isa(TriangulationStyle(trian),TriangulationStyle)
  bgtrian = get_background_triangulation(trian)
  @test isa(bgtrian,Triangulation)
  cell_id = get_cell_id(trian)
  @test isa(cell_id,AbstractArray) || isa(cell_id,SkeletonPair)
  cell_ref_map = get_cell_ref_map(trian)
  @test isa(cell_ref_map,AbstractArray) || isa(cell_ref_map,SkeletonPair)
end

# Some API

"""
    num_cells(trian::Triangulation) -> Int
"""
num_cells(trian::Triangulation) = length(get_cell_type(trian))

"""
    num_nodes(trian::Triangulation) -> Int
"""
num_nodes(trian::Triangulation) = length(get_node_coordinates(trian))

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
    is_first_order(trian::Triangulation) -> Bool
"""
function is_first_order(trian::Triangulation)
  reffes = get_reffes(trian)
  all(map(is_first_order,reffes))
end

"""
    get_cell_reffe(trian::Triangulation) -> Vector{<:LagrangianRefFE}

It is not desirable to iterate over the resulting array
for large number of cells if the underlying reference FEs
are of different Julia type.
"""
function get_cell_reffe(trian::Triangulation)
  type_to_reffe = get_reffes(trian)
  cell_to_type = get_cell_type(trian)
  expand_cell_data(type_to_reffe,cell_to_type)
end

"""
"""
function get_cell_ref_coordinates(trian::Triangulation)
  type_to_reffe = get_reffes(trian)
  type_to_coords = map(get_node_coordinates,type_to_reffe)
  cell_to_type = get_cell_type(trian)
  expand_cell_data(type_to_coords,cell_to_type)
end

"""
    get_cell_shapefuns(trian::Triangulation) -> Vector{<:Field}
"""
function get_cell_shapefuns(trian::Triangulation)
  type_to_reffes = get_reffes(trian)
  cell_to_type = get_cell_type(trian)
  type_to_shapefuns = map(get_shapefuns, type_to_reffes)
  expand_cell_data(type_to_shapefuns,cell_to_type)
end

"""
    get_cell_map(trian::Triangulation) -> Vector{<:Field}
"""
function get_cell_map(trian::Triangulation)
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_shapefuns = get_cell_shapefuns(trian)
  lazy_map(linear_combination,cell_to_coords,cell_to_shapefuns)
end

function get_node_coordinates(trian::Triangulation)
  @notimplemented
end

function get_cell_nodes(trian::Triangulation)
  @notimplemented
end

#"""
#"""
#function CellField(object,trian::Triangulation)
#  cm = get_cell_map(trian)
#  convert_to_cell_field(object,cm)
#end
#

# Helpers for Triangulation

#"""
#    restrict(cf::CellField,trian::Triangulation)
#"""
#function restrict(cf::CellField,trian::Triangulation)
#  _cf = to_ref_space(cf)
#  a = get_array(_cf)
#  r = restrict(a,trian)
#  axs = reindex(get_cell_axes(cf),trian)
#  _restrict_cell_field(r,axs,MetaSizeStyle(cf),trian)
#end
#
#function _restrict_cell_field(r::AbstractArray,axs::AbstractArray,msize_style::Val,trian)
#  cm = get_cell_map(trian)
#  GenericCellField(r,cm,Val(true),axs,msize_style)
#end
#
#function _restrict_cell_field(r::SkeletonPair,axs::SkeletonPair,msize_style::Val,trian)
#  cm = get_cell_map(trian)
#  la = r.plus
#  ra = r.minus
#  l = GenericCellField(la,cm,Val(true),axs.plus,msize_style)
#  r = GenericCellField(ra,cm,Val(true),axs.minus,msize_style)
#  merge_cell_fields_at_skeleton(l,r)
#end
#
#"""
#    CellQuadrature(trian::Triangulation, degree::Integer)
#"""
#function CellQuadrature(trian::Triangulation, degree::Integer)
#  polytopes = map(get_polytope,get_reffes(trian))
#  cell_type = get_cell_type(trian)
#  CellQuadrature(degree,polytopes,cell_type)
#end
#
#"""
#    integrate(cell_field,trian::Triangulation,quad::CellQuadrature)
#
#The `cell_field` is aligned with the cells in `trian`
#"""
#function integrate(cell_field,trian::Triangulation,quad::CellQuadrature)
#  cell_map = get_cell_map(trian)
#  integrate(cell_field,cell_map,quad)
#end
#
#function CellField(value::Number,trian::Triangulation,quad::CellQuadrature)
#  CellField(value,get_cell_map(trian),quad)
#end

