
"""
Abstract type that represents both a `CellField` and a `CellBasis`
"""
abstract type EvaluableCellArray{D,T,N} end

evaluate(::EvaluableCellArray{D,T,N} where {D,T,N},::CellPoints{D} where D)::CellArray{T,N} = @abstractmethod

# Concrete implementations

"""
Implements the results of a binary operation between two instances of `EvaluableCellArray`
"""
struct EvaluableCellArrayFromBinaryOp{D,T,N,A,B} <: EvaluableCellArray{D,T,N}
  a::EvaluableCellArray{D,T,A}
  b::EvaluableCellArray{D,T,B}
  op
end

function evaluate(self::EvaluableCellArrayFromBinaryOp,points::CellPoints)
  avals = evaluate(self.a,points)
  bvals = evaluate(self.b,points)
  op(avals,bvals)
end
