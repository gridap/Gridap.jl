module CellQuadratures

# Dependencies of this module

using Gridap
using Gridap.Helpers

# Functionality provided by this module

export CellQuadrature
export ConstantCellQuadrature
import Gridap: coordinates
import Gridap: weights

# Abstract types and interfaces

"""
Abstract type representing a collection of quadratures, one for each cell
"""
const CellQuadrature{D} = CellValue{<:Quadrature{D}}

coordinates(::CellQuadrature{D} where D )::CellPoints{D} = @abstractmethod

weights(::CellQuadrature)::CellValues{Float64} = @abstractmethod

function CellQuadrature end

# Pretty printing

import Base: show

function show(io::IO,self::CellQuadrature{D}) where {D}
  print(io,"CellQuadrature object")
end

function show(io::IO,::MIME"text/plain",quad::CellQuadrature{D}) where {D}
  show(io,quad)
  print(io,":")
  print(io,"\n dim: $D")
  print(io,"\n ncells: $(length(quad))")
end

# Concrete structs

"""
A concrete implementation of CellQuadrature for the particular case
that all cells have the same quadrature
"""
const ConstantCellQuadrature{D} = ConstantCellValue{Quadrature{D}}

function ConstantCellQuadrature(quad::Quadrature{D}, l::Int) where D
  ConstantCellQuadrature{D}(quad, l)
end

function coordinates(self::ConstantCellQuadrature)
  ConstantCellValue(coordinates(self.value), self.length)
end

function weights(self::ConstantCellQuadrature)
  ConstantCellValue(weights(self.value), self.length)
end

end # module CellQuadratures
