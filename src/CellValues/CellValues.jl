module CellValues

using Test
using Gridap
using Gridap.Helpers

export CellValue
export IterCellValue
export IndexCellValue

import Base: iterate
import Base: length
import Base: eltype
import Base: size
import Base: getindex
import Base: IndexStyle
import Base: show

export test_iter_cell_value
export test_index_cell_value

abstract type IterCellValue{T} end

function iterate(::IterCellValue{T})::Union{Nothing,Tuple{T,Any}} where T
  @abstractmethod
end

function iterate(::IterCellValue{T},state)::Union{Nothing,Tuple{T,Any}} where T
  @abstractmethod
end

length(::IterCellValue)::Int = @abstractmethod

eltype(::Type{C}) where C <: IterCellValue{T} where T = T

# Indexable cell Values

abstract type IndexCellValue{T,N} <: AbstractArray{T,N} end

function getindex(::IndexCellValue{T,N}, ::Vararg{<:Integer,N})::T where {T,N}
  @abstractmethod
end

size(x::IndexCellValue) = @abstractmethod

IndexStyle(::Type{<:IndexCellValue{T,N}} where {T,N}) = IndexLinear()

# Cell Values

const CellValue{T} = Union{IterCellValue{T},IndexCellValue{T}}

function show(io::IO,self::CellValue)
  for (i, a) in enumerate(self)
    println(io,"$i -> $a")
  end
end

function show(io::IO,::MIME"text/plain",self::CellValue)
  show(io,self)
end

# Testers

function test_iter_cell_value(icv::CellValue{T},a::AbstractArray{T}) where T
  _test_iter_cell_value(icv,a)
end

function _test_iter_cell_value(icv,a)

  @test length(icv) == length(a)

  i = 0
  for v in icv
    i += 1
    @assert _eq(v,a[i])
  end

  @test i == length(a)

  for v in icv
    @assert eltype(icv) == typeof(v)
  end

end

_eq(a,b) = a == b

_eq(a::Real,b::Real) = a ≈ b

function test_index_cell_value(icv::IndexCellValue{T},a::AbstractArray{T}) where T
  _test_index_cell_value(icv,a)
end

function _test_index_cell_value(icv,a)

  @test size(icv) == size(a)

  for i in eachindex(icv)
    v = icv[i]
    @assert v == a[i]
  end

  @test IndexStyle(icv) == IndexStyle(a)

end

end # module CellValues
