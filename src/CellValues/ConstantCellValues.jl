module ConstantCellValues

using Gridap

export ConstantCellValue
export ConstantCellNumber
export ConstantCellArray
export ConstantCellVector
export ConstantCellMatrix
export ConstantCellMap

import Gridap: apply
import Gridap: reindex
import Gridap: evaluate
import Base: size
import Base: getindex
import Base: broadcast
import Base: ==, ≈

struct ConstantCellValue{T} <: IndexCellValue{T,1}
  value::T
  length::Int
end

size(self::ConstantCellValue) = (self.length,)

function getindex(c::ConstantCellValue,i::Integer)
  c.value
end

const ConstantCellNumber = ConstantCellValue{<:NumberLike}

const ConstantCellArray{T,N} = ConstantCellValue{Array{T,N}}

const ConstantCellVector{T} = ConstantCellArray{T,1}

const ConstantCellMatrix{T} = ConstantCellArray{T,2}

const ConstantCellMap{S,M,T,N} = ConstantCellValue{<:Map{S,M,T,N}}

function ConstantCellArray(v::AbstractArray{T,N},l::Integer) where {T,N}
  ConstantCellArray{T,N}(v,l)
end

function ConstantCellVector(v::AbstractVector{T},l::Integer) where T
  ConstantCellVector{T}(v,l)
end

function ConstantCellMatrix(v::AbstractMatrix{T},l::Integer) where T
  ConstantCellMatrix{T}(v,l)
end

function ConstantCellMap(v::Map,l::Integer)
  ConstantCellValue(v,l)
end

function ConstantCellNumber(v::NumberLike,l::Integer)
  ConstantCellValue(v,l)
end

function (==)(a::ConstantCellNumber,b::ConstantCellNumber)
  _eq_kernel(==,a,b)
end

function (≈)(a::ConstantCellNumber,b::ConstantCellNumber)
  _eq_kernel(≈,a,b)
end

function (==)(a::ConstantCellArray,b::ConstantCellArray)
   _eq_kernel(==,a,b)
end

function (≈)(a::ConstantCellArray,b::ConstantCellArray)
  _eq_kernel(≈,a,b)
end

function _eq_kernel(op,a,b)
  !(op(a.value,b.value)) && return false
  length(a) != length(b) && return false
  return true
end

function apply(k::NumberKernel,v::Vararg{<:ConstantCellValue})
  vals = [vi.value for vi in v]
  l = _compute_l(v...)
  val = compute_value(k,vals...)
  ConstantCellNumber(val,l)
end

function apply(k::ArrayKernel,v::Vararg{<:ConstantCellValue})
  vals = [vi.value for vi in v]
  l = _compute_l(v...)
  val = compute_value(k,vals...)
  ConstantCellArray(val,l)
end

function apply(k::ArrayKernel,m::ConstantCellMap)
  w = apply(k,m.value)
  ConstantCellMap(w,m.length)
end

function apply(k::ArrayKernel,m::ConstantCellMap,v::Vararg{<:ConstantCellValue})
  vals = [vi.value for vi in v]
  l = _compute_l(m,v...)
  w = apply(k,m.value,vals...)
  ConstantCellMap(w,l)
end

function evaluate(cm::ConstantCellMap{S,M},ca::ConstantCellArray{<:S,M}) where {S,M}
  @assert length(cm) == length(ca)
  m = cm.value
  a = ca.value
  r = evaluate(m,a)
  ConstantCellArray(r,length(cm))
end

function _compute_l(v...)
  @assert length(v) > 0
  v1, = v
  l1 = length(v1)
  @assert all([ length(vi) == l1 for vi in v ])
  l1
end

function reindex(values::ConstantCellValue, indices::CellValue{<:IndexLike})
  ConstantCellValue(values.value,length(indices))
end

function reindex(values::ConstantCellValue, indices::IndexCellValue{<:IndexLike})
  ConstantCellValue(values.value,length(indices))
end

broadcast(f,cvs::ConstantCellValue) = ConstantCellValue(f(cvs.value),cvs.length)

end # module ConstantCellValues
