
"""
    abstract type CellMap <: GridapType end

Abstract type representing an object `ϕ` that maps a `CellPoint` object `q` to
a `CellPoint` object `x=ϕ(q)`.
It is represented as an array of `Field` objects, plus an indirection array of indices.
"""
abstract type CellMap <: GridapType end

"""
"""
function get_array(a::CellMap)
  @abstractmethod
end

"""
"""
function get_cell_id(a::CellMap)
  @abstractmethod
end

"""
Number of items in the range of `get_cell_id`
"""
function num_cell_ids(a::CellMap)
  @abstractmethod
end

function get_memo(a::CellMap)
  @abstractmethod
end

"""
"""
function test_cell_map(ϕ::CellMap)
  @test isa(get_array(ϕ),AbstractArray{<:Field})
  @test isa(get_cell_id(ϕ),AbstractArray{<:Integer})
  @test isa(num_cell_ids(ϕ),Int)
end

"""
"""
struct GenericCellMap <: CellMap
  array::AbstractArray{<:Field}
  cell_id::AbstractArray{<:Integer}
  num_cell_ids::Int
  memo::Dict
  function GenericCellMap(
    array::AbstractArray{<:Field}, cell_id::AbstractArray{<:Integer}, num_cell_ids::Int)
    new(array,cell_id,num_cell_ids,Dict())
  end
end

"""
"""
function GenericCellMap(array::AbstractArray{<:Field})
  n = length(array)
  cell_id = identity_vector(n)
  GenericCellMap(array,cell_id,n)
end

get_array(x::GenericCellMap) = x.array
get_cell_id(x::GenericCellMap) = x.cell_id
num_cell_ids(x::GenericCellMap) = x.num_cell_ids
get_memo(x::GenericCellMap) = x.memo

# Default API

Base.length(ϕ::CellMap) = length(get_array(ϕ))

"""
    evaluate(f::CellMap,x::CellPoint)
"""
function evaluate(f::CellMap,x::CellPoint)
  key = (:evaluate,objectid(x))
  memo = get_memo(f)
  if !haskey(memo,key)
    memo[key] = compute_evaluate(f,x)
  end
  memo[key]
end

function compute_evaluate(ϕ::CellMap,q::CellPoint)
  ϕids = get_cell_id(ϕ)
  qids = get_cell_id(q)
  array = evaluate_field_array(reindex(get_array(ϕ),qids),get_array(q))
  cell_ids = reindex(ϕids,qids)
  x = GenericCellPoint(array,cell_ids,num_cell_ids(ϕ))
  MappedCellPoint(x,ϕ,q)
end

"""
    (f::CellMap)(x)

Functor like-evaluation `f(x)` for `CellMap` objects. Syntactic sugar for `evaluate(f,x)`.
"""
function (f::CellMap)(x)
  evaluate(f,x)
end

struct MappedCellPoint <: CellPoint
  x::CellPoint
  ϕ::CellMap
  q::CellPoint
end

get_array(x::MappedCellPoint) = get_array(x.x)
get_cell_id(x::MappedCellPoint) = get_cell_id(x.x)
num_cell_ids(x::MappedCellPoint) = num_cell_ids(x.x)

"""
    Base.:(==)(a::CellMap,b::CellMap)

Check if two `CellMap` objects are the same. Uses `===` by default. 
"""
function Base.:(==)(a::CellMap,b::CellMap)
  if length(get_array(a)) != length(get_array(b))
    return false
  end
  a === b || get_array(a) === get_array(b) || get_array(a) == get_array(b)
end

"""
    gradient(cf::CellMap)
"""
function gradient(cf::CellMap)
  key = :gradient
  memo = get_memo(cf)
  if !haskey(memo,key)
    memo[key] = compute_gradient(cf)
  end
  memo[key]
end

function compute_gradient(ϕ::CellMap)
  cf = GenericCellField(get_array(ϕ))
  gradient(cf)
end

function Base.:∘(a::CellMap,b::CellMap)
  aids = get_cell_id(a)
  bids = get_cell_id(b)
  _a = reindex(get_array(a),bids)
  array = compose_field_arrays(_a,get_array(b))
  ids = reindex(aids,bids)
  GenericCellMap(array,ids,num_cell_ids(a))
end

"""
    inverse_map(ϕ::CellMap)

Return another `CellMap` representing the inverse map of the given map `ϕ`.
By defalut, it returns an instance of `InverseCellMap`, which keeps the direct cell map `ϕ` as metadata.
"""
inverse_map(ϕ::CellMap) = InverseCellMap(ϕ)

inverse_map(ϕ::AbstractArray) = @notimplemented

"""
    struct InverseCellMap <: CellField
      ϕ::CellField
    end

Struct representing the inverse of `ϕ` that keeps the direct cell map `ϕ` as metadata.
"""
struct InverseCellMap <: CellMap
  ϕ::CellMap
end

function get_array(ϕinv::InverseCellMap)
  ids = get_cell_id(ϕinv.ϕ)
  @notimplementedif ! isa(ids,IdentityVector)
  inverse_map(get_array(ϕinv.ϕ))
end

function get_cell_id(ϕinv::InverseCellMap)
  ids = get_cell_id(ϕinv.ϕ)
  @notimplementedif ! isa(ids,IdentityVector)
  ids
end

Base.length(ϕinv::InverseCellMap) = num_cell_ids(ϕinv.ϕ)
num_cell_ids(ϕinv::InverseCellMap) = length(ϕinv.ϕ)

#function Arrays.reindex(f::InverseCellMap,a::AbstractVector)
#  InverseCellMap(reindex(f.direct_cell_map,a),reindex(f.inverse_cell_map))
#end

"""
face_map = cell_map∘refface_to_refcell_map
"""
struct FaceMap <: CellMap
  face_map::CellMap
  cell_map::CellMap
  refface_to_refcell_map::CellMap
  memo::Dict
  function FaceMap(
  face_map::CellMap, cell_map::CellMap, refface_to_refcell_map::CellMap)
    new(face_map,cell_map,refface_to_refcell_map,Dict())
  end
end

get_array(a::FaceMap) = get_array(a.face_map)
get_cell_id(a::FaceMap) = get_cell_id(a.face_map)
num_cell_ids(a::FaceMap) = num_cell_ids(a.face_map)
get_memo(a::FaceMap) = a.memo





