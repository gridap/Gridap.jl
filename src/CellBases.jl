export CellBasis
export evaluate, gradient, mapderivatives

"""
Abstract type that represents a cell-wise basis for a field space,
where T is the type of value and D the dimension of the domain
"""
const CellBasis{D,T} = EvaluableCellArray{D,T,2} where {D,T}

evaluate(::CellBasis{D,T} where {D,T} ,::CellPoints{D} where D)::CellBasisValues{T}= @abstractmethod

"""
Returns another CellBasis object that represents the gradient
TG is a value whose rank is one order grater than the one of T
"""
gradient(::CellBasis{D,T} where {D,T})::CellBasis{D,TG} = @abstractmethod

"""
Returns another CellBasis, whose spatial
derivatives are respect to the coordinates of
the range space of geomap
"""
function mapderivatives(
  self::CellBasis{D,T}, geomap::CellField{D,Point{D}} ) where {D,T}
  CellBasisWithMappedDerivatives(self,geomap)
end
