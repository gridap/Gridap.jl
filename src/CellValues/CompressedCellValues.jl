module CompressedCellValues

using Gridap

export IterCompressedCellValue
export IndexCompressedCellValue

import Base: iterate
import Base: length
import Base: size
import Base: getindex

struct IterCompressedCellValue{T,A} <: IterCellValue{T}
  values::Vector{T}
  ptrs::A

  function IterCompressedCellValue(
    values::Vector{T}, ptrs) where T
    @assert hasmethod(iterate,Tuple{typeof(ptrs)})
    @assert hasmethod(length,Tuple{typeof(ptrs)})
    A = typeof(ptrs)
    new{T,A}(values,ptrs)
  end

end

@inline function iterate(cv::IterCompressedCellValue)
  inext = iterate(cv.ptrs)
  _iterate(cv,inext)
end

@inline function iterate(cv::IterCompressedCellValue,state)
  inext = iterate(cv.ptrs,state)
  _iterate(cv,inext)
end

@inline function _iterate(cv,inext)
  if inext === nothing; return nothing; end
  i, istate = inext
  v = cv.values[i]
  (v,istate)
end

length(cv::IterCompressedCellValue) = length(cv.ptrs)

struct IndexCompressedCellValue{T,A} <: IndexCellValue{T,1}
  values::Vector{T}
  ptrs::A

  function IndexCompressedCellValue(
    values::Vector{T}, ptrs::AbstractArray) where T
    A = typeof(ptrs)
    new{T,A}(values,ptrs)
  end

end

function IndexCompressedCellValue(cv::ConstantCellValue)
  values = [cv.value,]
  ptrs = ConstantCellValue(1,length(cv))
  IndexCompressedCellValue(values,ptrs)
end

function getindex(
  cv::IndexCompressedCellValue{T,A}, i::Integer) where {T,A}
  j = cv.ptrs[i]
  cv.values[j]
end

size(cv::IndexCompressedCellValue) = (length(cv.ptrs),)

end # module
