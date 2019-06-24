module FindLocalIndex

using Gridap

import Base: size
import Base: getindex

export find_local_index
export get_local_item

function find_local_index(gids::CellNumber,lids_to_gids::CellArray)
  LocalIndices(Int,gids,lids_to_gids)
end

function get_local_item(a::IndexCellArray,litem::Int)
  LocalItems(a,litem)
end

struct LocalIndices{T,A,B} <: IndexCellNumber{T,1}
  gids::A
  lids_to_gids::B
end

function LocalIndices(
  T::Type,gids::IndexCellNumber,lids_to_gids::IndexCellArray)
  @assert T <: Integer
  A = typeof(gids)
  B = typeof(lids_to_gids)
  LocalIndices{T,A,B}(gids,lids_to_gids)
end

size(l::LocalIndices) = (length(l.gids),)

function getindex(l::LocalIndices{T},i::Integer) where T
  j = l.gids[i]
  lid_to_gid = l.lids_to_gids[j]
  k = zero(T)
  for g in lid_to_gid
    k += 1
    if g == i
      return k
    end
  end
  return zero(T)
end

struct LocalItems{T,A} <: IndexCellNumber{T,1}
  a::A
  lindex::Int
end

function LocalItems(a::IndexCellArray{T},lindex::Int) where T
  A = typeof(a)
  LocalItems{T,A}(a,lindex)
end

size(l::LocalItems) = (length(l.a),)

function getindex(l::LocalItems,i::Integer)
  v = l.a[i]
  v[l.lindex]
end

end # module
