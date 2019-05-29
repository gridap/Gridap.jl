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

