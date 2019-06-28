module CompressedCellValues

using Gridap

export IterCompressedCellValue

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

end # module
