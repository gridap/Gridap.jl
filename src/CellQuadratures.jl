module CellQuadratures

# Dependencies of this module

using Numa.Helpers
using Numa.FieldValues
using Numa.Quadratures
using Numa.CellValues
using Numa.CellMaps
using Numa.Geometry

# Functionality provided by this module

export CellQuadrature
export ConstantCellQuadrature
import Base.Iterators: zip
import Numa: coordinates
import Numa: weights
import Numa: quadrature

# Abstract types and interfaces

# @santiagobadia : I would create a JoinCellValues struc, and define this one as
# sub-type with nothing but coordinates and weights methods... Do we really need
# it?

"""
Abstract type representing a collection of quadratures, one for each cell
"""
abstract type CellQuadrature{D} end

coordinates(::CellQuadrature{D} where D )::CellPoints{D} = @abstractmethod

weights(::CellQuadrature)::CellValues{Float64} = @abstractmethod

function zip(self::CellQuadrature)
  c = coordinates(self)
  w = weights(self)
  zip(c,w)
end

# Factories

"""
Factory function to create CellQuadrature objects in a convenient way
"""
function quadrature(trian::Triangulation;order::Int)
  _quadrature(celltypes(trian),order)
end

# Concrete structs

"""
A concrete implementation of CellQuadrature for the particular case
that all cells have the same quadrature
"""
struct ConstantCellQuadrature{D} <: CellQuadrature{D}
  coords::ConstantCellArray{Point{D},1}
  weights::ConstantCellArray{Float64,1}
end

function ConstantCellQuadrature(c::Array{Point{D},1},w::Array{Float64,1},l::Int) where D
  @assert length(c) == length(w)
  coords = ConstantCellValue(c,l)
  weights = ConstantCellValue(w,l)
  ConstantCellQuadrature{D}(coords,weights)
  # santiagobadia : Be careful here... without D it does not work because
  # ConstantCellValue not templatized by dim. Why did it work with
  # ConstantCellArray
end

function ConstantCellQuadrature(quad::Quadrature{D} where D,l::Int)
  c = coordinates(quad)
  w = weights(quad)
  ConstantCellQuadrature(c,w,l)
end

coordinates(self::ConstantCellQuadrature) = self.coords

weights(self::ConstantCellQuadrature) = self.weights

# Helpers

_quadrature(ct,order) = @notimplemented

function _quadrature(ct::ConstantCellValue{NTuple{Z,Int}},order) where Z
  t = celldata(ct)
  q = quadrature(t,order=order)
  ConstantCellQuadrature(q,length(ct))
end

end # module CellQuadratures
