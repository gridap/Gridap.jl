
"""
Abstract type that represents both a `CellField` and a `CellBasis`
"""
abstract type EvaluableCellArray{D,T,N} end

evaluate(::EvaluableCellArray{D,T,N} where {D,T,N},::CellPoints{D} where D)::CellArray{T,N} = @abstractmethod

# Concrete implementations

"""
Implements the results of a binary operation `op` between two instances `a` and `b` of `EvaluableCellArray`
"""
struct EvaluableCellArrayFromBinaryOp{D,T,N,A,B} <: EvaluableCellArray{D,T,N}
  a::EvaluableCellArray{D,T,A}
  b::EvaluableCellArray{D,T,B}
  op
end

"""
Compute the operation `op(a,b)` for a `EvaluableCellArrayFromBinaryOp`
"""
function evaluate(self::EvaluableCellArrayFromBinaryOp,points::CellPoints)
  avals = evaluate(self.a,points)
  bvals = evaluate(self.b,points)
  self.op(avals,bvals)
end

 """
 Implements the composition `a ∘ b` of two instances `a` and `b` of `EvaluableCellArray`
 """
 struct EvaluableCellArrayFromComposition{D,T,N} <: EvaluableCellArray{D,T,N}
   a::EvaluableCellArray{D,T,N}
   b::EvaluableCellArray{D,Point{D},1}
 end
# @santiagobadia : Not sure this method will be pf practical use in our FE integration
# machinery, since we cannot provide the more concrete type. Should we create a composition
# for every concrete CellBasis and CellField? Or simply checking that D,T,N are right...

 Base.:∘(f::EvaluableCellArray{D,T,N},
         g::EvaluableCellArray{D,Point{D},1}) where{D,T,N} =
		 EvaluableCellArrayFromComposition(f,g)

"""
Evaluate a `EvaluableCellArrayFromComposition`
"""
function evaluate(self::EvaluableCellArrayFromComposition, points::CellPoints)
  gpoins = evaluate(self.g,points)
  bvals = evaluate(self.b,gpoins)
end
