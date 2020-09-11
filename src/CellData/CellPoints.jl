
"""
    abstract type CellPoint <: GridapType end

Abstract type representing evaluation points in the cells (faces) of a mesh.
It is a plain array of array of points plus an indirection array of indices, which
allows to evaluate only a subset of cells.
"""
abstract type CellPoint <: GridapType end

"""
"""
function get_array(a::CellPoint)
  @abstractmethod
end

"""
"""
function get_cell_id(a::CellPoint)
  @abstractmethod
end

"""
Number of items in the range of `get_cell_id`
"""
function num_cell_ids(a::CellPoint)
  @abstractmethod
end

function test_cell_point(x::CellPoint)
  @test isa(get_array(x),AbstractArray{<:AbstractArray{<:Point}})
  @test isa(get_cell_id(x),AbstractArray{<:Integer})
  @test isa(num_cell_ids(x),Int)
end

Base.length(x::CellPoint) = length(get_array(x))

#function Arrays.reindex(x::CellPoint,a::AbstractArray)
#  @notimplementedif !isa(get_cell_id(x),IdentityVector)
#  array = reindex(get_array(x),a)
#  GenericCellPoint(array,a,length(x))
#end

struct GenericCellPoint <: CellPoint
  array::AbstractArray{<:AbstractArray{<:Point}}
  cell_id::AbstractArray{<:Integer}
  num_cell_ids::Int
end

function GenericCellPoint(array::AbstractArray{<:AbstractArray{<:Point}})
  n = length(array)
  cell_id = identity_vector(n)
  GenericCellPoint(array,cell_id,n)
end

get_array(x::GenericCellPoint) = x.array
get_cell_id(x::GenericCellPoint) = x.cell_id
num_cell_ids(x::GenericCellPoint) = x.num_cell_ids

# Lazy append

function Arrays.lazy_append(a::CellPoint,b::CellPoint)
  @notimplementedif !isa(get_cell_id(a),IdentityVector)
  @notimplementedif !isa(get_cell_id(b),IdentityVector)
  array = lazy_append(get_array(a),get_array(b))
  GenericCellPoint(array)
end

