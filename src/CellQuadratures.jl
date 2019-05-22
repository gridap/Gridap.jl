module CellQuadratures

# Dependencies of this module

using Gridap.Helpers
using Gridap.FieldValues
using Gridap.Quadratures
using Gridap.CellValues
using Gridap.CellValues.ConstantCellValues
using Gridap.CellMaps
using Gridap.Geometry

# Functionality provided by this module

export CellQuadrature
export ConstantCellQuadrature
import Base.Iterators: zip
import Gridap: coordinates
import Gridap: weights

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
function CellQuadrature(trian::Triangulation;order::Int)
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
  q = Quadrature(t,order=order)
  ConstantCellValue(q,length(ct))
  ConstantCellQuadrature(q,length(ct))
end

end # module CellQuadratures
