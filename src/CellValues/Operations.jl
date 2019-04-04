
# Unary operations on CellValue

for op in (:+, :-, :(inv), :(det))
  @eval begin
    function ($op)(a::CellValue)
      CellValueFromUnaryOp($op,a)
    end
  end
end

struct CellValueFromUnaryOp{T,O<:Function,C<:CellValue} <: IterCellValue{T}
  op::O
  values::C
end

function CellValueFromUnaryOp(op::Function,values::CellValue{T}) where T
  O = typeof(op)
  C = typeof(values)
  S = Base._return_type(op,Tuple{T})
  CellValueFromUnaryOp{S,O,C}(op,values)
end

length(self::CellValueFromUnaryOp) = length(self.values)

@inline function iterate(self::CellValueFromUnaryOp)
  next = iterate(self.values)
  if next === nothing; return nothing end
  a, state = next
  (self.op(a), state)
end

@inline function iterate(self::CellValueFromUnaryOp,state)
  next = iterate(self.values,state)
  if next === nothing; return nothing end
  a, state = next
  (self.op(a), state)
end

#struct CellValueFromBinOp{O<:Function,A<:CellValue,B<:CellValue,T} <: CellValue{T}
#  op::O
#  a::A
#  b::B
#end
#
#function CellValueFromBinOp(op::Function,a::CellValue{T},b::CellValue{S}) where {T,S}
#  O = typeof(op)
#  A = typeof(a)
#  B = typeof(b)
#  R = Base._result_type(op,Tuple{T,S})
#  CellValueFromBinOp{O,C,R}(op,values)
#end
#
#@inline function iterate(self::CellValueFromUnaryOp)
#  anext = iterate(self.values)
#  if anext === nothing; return nothing end
#  if bnext === nothing; return nothing end
#  a, astate = anext
#  b, bstate = bnext
#  state = (astate,bstate)
#  (self.op(a,b), state)
#end
#
#@inline function iterate(self::CellValueFromUnaryOp,state)
#  astate, bstate = state
#  anext = iterate(self.a,astate)
#  if anext === nothing; return nothing end
#  a, astate = anext
#  (self.op(a), state)
#end


