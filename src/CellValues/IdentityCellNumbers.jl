module IdentityCellNumbers

using Gridap

export IdentityCellNumber
import Gridap: reindex
import Base: size
import Base: getindex

struct IdentityCellNumber{T} <: IndexCellValue{T,1}
  length::Int
end

function IdentityCellNumber(::Type{T},l::Integer) where T <: Integer
  IdentityCellNumber{T}(l)
end

function getindex(c::IdentityCellNumber{T},i::Integer) where T
  @assert i > 0
  @assert i <= c.length
  j::T = i
  j
end

size(c::IdentityCellNumber) = (c.length,)

function reindex(values::IndexCellValue, indices::IdentityCellNumber{T}) where T
  @assert length(values) == length(indices)
  values
end

function reindex(values::IndexCellNumber, indices::IdentityCellNumber{T}) where T
  @assert length(values) == length(indices)
  values
end

function reindex(values::IndexCellArray, indices::IdentityCellNumber{T}) where T
  @assert length(values) == length(indices)
  values
end

end # module
