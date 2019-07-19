module NonIterableCellMaps

using Gridap
using Gridap.Helpers

export UnimplementedMap
export NonIterableCellMap
import Base: iterate
import Gridap: evaluate
import Gridap: return_size

struct UnimplementedMap{S,M,T,N} <: Map{S,M,T,N} end

function evaluate!(
  this::UnimplementedMap{S,M,T,N},
  points::AbstractArray{<:S,M},
  v::AbstractArray{T,N}) where {S,M,T,N}
  @notimplemented
end

function return_size(::UnimplementedMap{S,M,T,N},::NTuple{M,Int}) where {S,M,T,N}
  @notimplemented
end

abstract type NonIterableCellMap{S,M,T,N} <: IterCellValue{UnimplementedMap{S,M,T,N}} end

function iterate(::NonIterableCellMap)
  _error()
end

function iterate(::NonIterableCellMap,state)
  _error()
end

function _error()
  error("Iteration is intentionally disabled for this type.")
end

end # module
