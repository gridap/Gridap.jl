module CellQuadratures

# Dependencies of this module

using Numa.Helpers
using Numa.FieldValues
using Numa.Quadratures
using Numa.CellValues
using Numa.CellValues.ConstantCellValues
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

"""
Abstract type representing a collection of quadratures, one for each cell
"""
const CellQuadrature{D} = CellValue{Quadrature{D}}

coordinates(::CellQuadrature{D} where D )::CellPoints{D} = @abstractmethod

weights(::CellQuadrature)::CellValues{Float64} = @abstractmethod

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
const ConstantCellQuadrature{D} = ConstantCellValue{Quadrature{D}}

ConstantCellQuadrature(quad::Quadrature{D}, l::Int) where D =  ConstantCellQuadrature{D}(quad, l)

coordinates(self::ConstantCellQuadrature) = ConstantCellValue(coordinates(self.value), self.length)

weights(self::ConstantCellQuadrature) = ConstantCellValue(weights(self.value), self.length)

# Helpers

_quadrature(ct,order) = @notimplemented

function _quadrature(ct::ConstantCellValue{NTuple{Z,Int}},order) where Z
  t = celldata(ct)
  q = quadrature(t,order=order)
  ConstantCellValue(q,length(ct))
  ConstantCellQuadrature(q,length(ct))
end

end # module CellQuadratures
