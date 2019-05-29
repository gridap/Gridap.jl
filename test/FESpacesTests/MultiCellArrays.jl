module MultiCellArrays

using Gridap
using Gridap.CachedArrays
using Gridap.CellValues

export MultiCellArray
export eachblock
import Base: iterate
import Base: getindex
import Base: length

struct MultiCellArray{T,N}
  cellarrays::Vector{<:CellValue{CachedArray{T,N,Array{T,N}}}}
  fieldids::Vector{NTuple{N,Int}}
end

# TODO create a constructor that checks that all return CachedArrays
# tread constant arrays efficiently (i.e., create ConstandCachedArray)

#TODO check at construction that we input at least one
# and that all have the same size
length(mca::MultiCellArray) = length(mca.cellarrays[1])

@inline function iterate(mca::MultiCellArray)
  zipped = zip(mca.cellarrays...)
  znext = iterate(zipped)
  if znext === nothing; return nothing end
  arrays, zstate = znext
  ma = MultiCachedArray(arrays,mca.fieldids)
  state = (ma, zstate, zipped)
  (ma, state)
end

@inline function iterate(mca::MultiCellArray,state)
  ma, zstate, zipped = state
  znext = iterate(zipped, zstate)
  if znext === nothing; return nothing end
  arrays, zstate = znext
  ma.arrays = arrays
  state = (ma, zstate, zipped)
  (ma, state)
end

mutable struct MultiCachedArray{L,T,N}
  arrays::NTuple{L,CachedArray{T,N,Array{T,N}}}
  fieldids::Vector{NTuple{N,Int}}
end

#Check that the length of arrays and fieldids is the same
length(ma::MultiCachedArray{L}) where L = L

getindex(ma::MultiCachedArray,i::Integer) = ma.arrays[i]

# TODO check if this allocates a new object
eachblock(ma::MultiCachedArray) = zip(ma.arrays,ma.fieldids)

end # module MultiCellArrays

