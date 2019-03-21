export CellField
export evaluate

abstract type CellField{D,T} end

evaluate(::CellField{D,T} where {D,T},::CellPoints{D} where D)::CellFieldValues{T} = @abstractmethod

"""
Returns another CellField object that represents the gradient.
TG is a value whose rank is one order grater than the one of T
"""
gradient(::CellField{D,T} where {D,T})::CellField{D,TG} = @abstractmethod

