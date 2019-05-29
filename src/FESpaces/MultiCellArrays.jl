module MultiCellArrays

using Gridap
using Gridap.CachedArrays
using Gridap.CellValues
using Gridap.CellValues.Operations: CellArrayFromBroadcastUnaryOp
using Gridap.CellValues.ConstantCellValues

export MultiCellArray
export eachblock
import Base: iterate
import Base: getindex
import Base: length

struct MultiCellArray{T,N}
  cellarrays::Vector{<:CellValue{CachedArray{T,N,Array{T,N}}}}
  fieldids::Vector{NTuple{N,Int}}
end

function MultiCellArray(cellarrays::Vector{<:CellArray{T,N}},fieldids::Vector{NTuple{N,Int}}) where {T,N}
  @assert length(cellarrays) > 0
  @assert all( [ length(ca) == length(cellarrays[1]) for ca in cellarrays ])
  _cellarrays = [ _prepare(ca) for ca in cellarrays]
  MultiCellArray{T,N}(_cellarrays,fieldids)
end

function _prepare(ca::CellValue{CachedArray{T,N,Array{T,N}}}) where {T,N}
  ca
end

function _prepare(ca::CellArray{T,N}) where {T,N}
  CellArrayFromBroadcastUnaryOp(+,ca)
end

function _prepare(ca::ConstantCellArray{T,N}) where {T,N}
  _ca = CachedArray(celldata(ca))
  ConstantCellValue(_ca,length(ca))
end

length(mca::MultiCellArray) = length(mca.cellarrays[1])

@inline function iterate(mca::MultiCellArray)
  zipped = zip(mca.cellarrays...)
  znext = iterate(zipped)
  if znext === nothing; return nothing end
  arrays, zstate = znext
  state = (zstate, zipped)
  (arrays, state)
end

@inline function iterate(mca::MultiCellArray,state)
  zstate, zipped = state
  znext = iterate(zipped, zstate)
  if znext === nothing; return nothing end
  arrays, zstate = znext
  state = (zstate, zipped)
  (arrays, state)
end

end # module MultiCellArrays

