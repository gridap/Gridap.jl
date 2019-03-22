export CellField
export evaluate

"""
Abstract type that represents a cell-wise field, where
T stands for the type that represents the field at a point
(e.g., scalar, vector, tensor) and D stands for the space
dimension
"""
const CellField{D,T} = EvaluableCellArray{D,T,1} where {D,T}

evaluate(::CellField{D,T} where {D,T},::CellPoints{D} where D)::CellFieldValues{T} = @abstractmethod

"""
Returns another CellField object that represents the gradient.
TG has a rank one order greater than the one of T
"""
gradient(::CellField{D,T} where {D,T})::CellField{D,TG} = @abstractmethod
