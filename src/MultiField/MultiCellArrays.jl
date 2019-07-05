module MultiCellArrays

using Gridap
using Gridap.CachedArrays
using Gridap.CellArrayApply: CellArrayFromKernel

export MultiCellArray
export MultiCellMatrix
export MultiCellVector
import Base: iterate
import Base: length
import Gridap: reindex

struct MultiCellArray{T,N}
  blocks::Vector{<:CellValue{CachedArray{T,N,Array{T,N}}}}
  fieldids::Vector{NTuple{N,Int}}
end

const MultiCellMatrix{T} = MultiCellArray{T,2}

const MultiCellVector{T} = MultiCellArray{T,1}

function MultiCellArray(cellarrays::Vector{<:CellArray{T,N}},fieldids::Vector{NTuple{N,Int}}) where {T,N}
  @assert length(cellarrays) > 0
  @assert all( [ length(ca) == length(cellarrays[1]) for ca in cellarrays ])
  _cellarrays = [ _prepare(ca) for ca in cellarrays]
  MultiCellArray{T,N}(_cellarrays,fieldids)
end

function MultiCellMatrix(cellmatrices::Vector{<:CellMatrix{T}},fieldids::Vector{NTuple{2,Int}}) where T
  MultiCellArray(cellmatrices,fieldids)
end

function MultiCellVector(cellvectors::Vector{<:CellVector{T}},fieldids::Vector{NTuple{1,Int}}) where T
  MultiCellArray(cellvectors,fieldids)
end

function _prepare(ca::CellValue{CachedArray{T,N,Array{T,N}}}) where {T,N}
  ca
end

function _prepare(ca::CellArray{T,N}) where {T,N}
  k = ArrayKernelFromBroadcastedFunction(+)
  CellArrayFromKernel(k,ca)
end

function _prepare(ca::ConstantCellArray{T,N}) where {T,N}
  _ca = CachedArray(ca.value)
  ConstantCellValue(_ca,length(ca))
end

length(mca::MultiCellArray) = length(mca.blocks[1])

@inline function iterate(mca::MultiCellArray)
  zipped = zip(mca.blocks...)
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

function reindex(mca::MultiCellArray,cn::CellNumber)
  blocks = [ reindex(block,cn)  for block in mca.blocks ]
  MultiCellArray(blocks,mca.fieldids)
end

end # module
