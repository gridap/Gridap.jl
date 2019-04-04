
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

# Binary operations on CellValue

for op in (:+, :-, :*, :/, :(inner), :(outer))
  @eval begin
    function ($op)(a::CellValue,b::CellValue)
      CellValueFromBinaryOp($op,a,b)
    end
  end
end

struct CellValueFromBinaryOp{T,O<:Function,A<:CellValue,B<:CellValue} <: IterCellValue{T}
  op::O
  a::A
  b::B
end

function CellValueFromBinaryOp(op::Function,a::CellValue{T},b::CellValue{S}) where {T,S}
  @assert length(a) == length(b)
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{T,S})
  CellValueFromBinaryOp{R,O,A,B}(op,a,b)
end

function length(self::CellValueFromBinaryOp)
  @assert length(self.a) == length(self.b)
  length(self.a)
end

@inline function iterate(self::CellValueFromBinaryOp)
  anext = iterate(self.a)
  bnext = iterate(self.b)
  if anext === nothing; return nothing end
  if bnext === nothing; return nothing end
  a, astate = anext
  b, bstate = bnext
  state = (astate,bstate)
  (self.op(a,b), state)
end

@inline function iterate(self::CellValueFromBinaryOp,state)
  astate, bstate = state
  anext = iterate(self.a,astate)
  bnext = iterate(self.b,bstate)
  if anext === nothing; return nothing end
  if bnext === nothing; return nothing end
  a, astate = anext
  b, bstate = bnext
  state = (astate,bstate)
  (self.op(a,b), state)
end

