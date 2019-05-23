module ConstantCellValues

using Gridap.CellValues
using Gridap.CellValues.Operations: cellsumsize

export ConstantCellValue
export ConstantCellArray
export ConstantCellVector
export ConstantCellMatrix

export celldata

import Gridap.CellValues: cellsum
import Gridap.CellValues: cellnewaxis
import Gridap.CellValues: cellmean #TODO
import Gridap: apply
import Base: +, -, *, /
import Base: ==
import LinearAlgebra: inv, det
import Gridap.FieldValues: inner, outer, meas

import Base: size
import Base: getindex
import Gridap.CellValues: cellsize

struct ConstantCellValue{T} <: IndexCellValue{T,1}
  value::T
  length::Int
end

celldata(self::ConstantCellValue) = self.value

size(self::ConstantCellValue) = (self.length,)

const ConstantCellArray{T,N} = ConstantCellValue{Array{T,N}}

const ConstantCellVector{T} = ConstantCellArray{T,1}

const ConstantCellMatrix{T} = ConstantCellArray{T,2}

function ConstantCellArray(v::AbstractArray{T,N},l::Integer) where {T,N}
  ConstantCellArray{T,N}(v,l)
end

function ConstantCellVector(v::AbstractVector{T},l::Integer) where T
  ConstantCellVector{T}(v,l)
end

function ConstantCellMatrix(v::AbstractMatrix{T},l::Integer) where T
  ConstantCellMatrix{T}(v,l)
end

cellsize(self::ConstantCellArray) = size(self.value)

function cellsum(self::ConstantCellArray{T,N};dim::Int) where {T,N}
  b = sum(self.value,dims=dim)
  s = cellsumsize(size(b),Val(dim))
  c = copy(reshape(b,s))
  ConstantCellValue(c,self.length)
end

function cellsum(self::ConstantCellArray{T,1};dim::Int) where T
  b = sum(self.value)
  ConstantCellValue(b,self.length)
end

function cellnewaxis(self::ConstantCellArray;dim::Int)
  s = [ v for v in size(self.value)]
  insert!(s,dim,1)
  shape = tuple(s...)
  c = copy(reshape(self.value,shape))
  ConstantCellValue(c,self.length)
end

getindex(self::ConstantCellValue,cell::Int) = celldata(self)

function (==)(a::ConstantCellValue,b::ConstantCellValue)
  celldata(a) != celldata(b) && return false
  length(a) != length(b) && return false
  return true
end

for op in (:+,:-,:*,:/,:(outer),:(inner))

  @eval begin
    function ($op)(a::ConstantCellValue,b::ConstantCellValue)
      @assert length(a) == length(b)
      c = _bin_op_kernel($op,celldata(a),celldata(b))
      ConstantCellValue(c,length(a))
    end
  end

end

function _bin_op_kernel(op,a,b)
  op(a,b)
end

function _bin_op_kernel(op,a::AbstractArray,b::AbstractArray)
  broadcast(op,a,b)
end

function _bin_op_kernel(op,a,b::AbstractArray)
  broadcast(op,a,b)
end

function _bin_op_kernel(op,a::AbstractArray,b)
  broadcast(op,a,b)
end

for op in (:+,:-,:(det),:(inv),:(meas))

  @eval begin
    function ($op)(a::ConstantCellValue)
      c = _unary_op_kernel($op,celldata(a))
      ConstantCellValue(c,length(a))
    end
  end

end

function _unary_op_kernel(op,a)
  op(a)
end

function _unary_op_kernel(op,a::AbstractArray)
  broadcast(op,a)
end

function apply(op::Function,a::ConstantCellValue)
  c = _unary_op_kernel(op,celldata(a))
  ConstantCellValue(c,length(a))
end

function apply(op::Function,a::ConstantCellValue,b::ConstantCellValue)
  @assert length(a) == length(b)
  c = _bin_op_kernel(op,celldata(a),celldata(b))
  ConstantCellValue(c,length(a))
end

end # module ConstantCellValues
