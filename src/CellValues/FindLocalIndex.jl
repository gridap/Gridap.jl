module FindLocalIndex

using Gridap

import Base: size
import Base: getindex

export find_local_index

function find_local_index(gids::CellNumber,lids_to_gids::CellArray)
  LocalIndices(Int,gids,lids_to_gids)
end

struct LocalIndices{T,A,B} <: IndexCellNumber{T,1}
  gids::A
  lids_to_gids::B
end

function LocalIndices(
  T::Type,gids::IndexCellNumber,lids_to_gids::IndexCellArray)
  @assert T <: Integer
  @assert length(gids) == length(lids_to_gids)
  A = typeof(gids)
  B = typeof(lids_to_gids)
  LocalIndices{T,A,B}(gids,lids_to_gids)
end

size(l::LocalIndices) = (length(l.gids),)

function getindex(l::LocalIndices{T},i::Integer) where T
  gid = l.gids[i]
  lid_to_gid = l.lids_to_gids[i]
  k = zero(T)
  for g in lid_to_gid
    k += 1
    if g == gid
      return k
    end
  end
  return zero(T)
end

end # module
