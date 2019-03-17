export Point, CellQuadrature, ConstantCellQuadrature
export coordinates, weights

using StaticArrays

"""
Type representing a point of D dimensions
"""
const Point{D} = SVector{D,Float64} where D

"""
Abstract type representing a collection of quadratures, one for each cell
"""
abstract type CellQuadrature{D} end

coordinates(::CellQuadrature{D} where D )::CellFieldValues{Point{D}} = @abstractmethod

weights(::CellQuadrature)::CellFieldValues{Float64} = @abstractmethod

function Base.Iterators.zip(self::CellQuadrature)
  c = coordinates(self)
  w = weights(self)
  zip(c,w)
end


# Concrete implementations

"""
A concrete implementation of CellQuadrature for the particular case
that all cells have the same quadrature
"""
struct ConstantCellQuadrature{D} <: CellQuadrature{D}
  coords::ConstantCellArray{Point{D},1}
  weights::ConstantCellArray{Float64,1}
end

function ConstantCellQuadrature(c::Array{Point{D},1} where D,w::Array{Float64,1},l::Int)
  @assert length(c) == length(w)
  coords = ConstantCellArray(c,l)
  weights = ConstantCellArray(w,l)
  ConstantCellQuadrature(coords,weights)
end

coordinates(self::ConstantCellQuadrature) = self.coords

weights(self::ConstantCellQuadrature) = self.weights
