
function Base.:+(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  CellArrayFromSum{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:-(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  CellArrayFromSub{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:*(a::CellArray{A,N},b::CellArray{B,N}) where {A,B,N}
  T = A * B
  CellArrayFromMul{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:/(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  CellArrayFromDiv{typeof(a),typeof(b),T,N}(a,b)
end

"""
Assumes that outer is defined by instances of T,S, and T,S as types as well
"""
function bouter(a::CellArray{T,N},b::CellArray{S,N}) where {T,S,N}
  R = outer(T,S)
  CellArrayFromOuter{typeof(a),typeof(b),R,N}(a,b)
end

"""
Assumes that inner is defined by instances of T and T as type as well
"""
function binner(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  R = inner(T,T)
  CellArrayFromInner{typeof(a),typeof(b),R,N}(a,b)
end

function Base.:(==)(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  length(a) != length(b) && return false
  cellsize(a) != cellsize(b) && return false
  for (ai,bi) in zip(a,b)
    ai != bi && return false
  end
  return true
end

"""
Assumes that det is defined for instances of T
and that the result is Float64
"""
function LinearAlgebra.det(self::CellArray{T,N}) where {T,N}
  CellArrayFromGivenUnaryOp{typeof(det),typeof(self),Float64,N}(self,det)
end

"""
Assumes that inv is defined for instances of T
"""
function LinearAlgebra.inv(self::CellArray{T,N}) where {T,N}
  CellArrayFromGivenUnaryOp{typeof(inv),typeof(self),T,N}(self,inv)
end

function cellsum(self::CellArray{T,N};dim::Int) where {T,N}
  CellArrayFromCellSum{dim,N-1,typeof(self),T}(self)
end

function cellnewaxis(self::CellArray{T,N};dim::Int) where {T,N}
  CellArrayFromCellNewAxis{dim,typeof(self),T,N+1}(self)
end

# Ancillary types associated with the operations above

struct CellArrayFromGivenUnaryOp{O<:Function,C,T,N} <: CellArrayFromElemUnaryOp{C,T,N}
  a::C
  op::O
end

inputcellarray(self::CellArrayFromGivenUnaryOp) = self.a

function computevals!(self::CellArrayFromGivenUnaryOp, a, v)
  v .= self.op.(a)
end

"""
Lazy sum of two cell arrays
"""
struct CellArrayFromSum{A,B,T,N} <: CellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::CellArrayFromSum) = self.a

rightcellarray(self::CellArrayFromSum) = self.b

function computevals!(::CellArrayFromSum, a, b, v)
  v .= a .+ b
end

"""
Lazy subtraction of two cell arrays
"""
struct CellArrayFromSub{A,B,T,N} <: CellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::CellArrayFromSub) = self.a

rightcellarray(self::CellArrayFromSub) = self.b

function computevals!(::CellArrayFromSub, a, b, v)
  v .= a .- b
end

"""
Lazy multiplication of two cell arrays
"""
struct CellArrayFromMul{A,B,T,N} <: CellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::CellArrayFromMul) = self.a

rightcellarray(self::CellArrayFromMul) = self.b

function computevals!(::CellArrayFromMul, a, b, v)
  v .= a .* b
end

"""
Lazy division of two cell arrays
"""
struct CellArrayFromDiv{A,B,T,N} <: CellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::CellArrayFromDiv) = self.a

rightcellarray(self::CellArrayFromDiv) = self.b

function computevals!(::CellArrayFromDiv, a, b, v)
  v .= a ./ b
end

"""
Lazy outer of two cell arrays
"""
struct CellArrayFromOuter{A,B,T,N} <: CellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::CellArrayFromOuter) = self.a

rightcellarray(self::CellArrayFromOuter) = self.b

function computevals!(::CellArrayFromOuter, a, b, v)
  v .= outer.(a,b)
end

"""
Lazy inner of two cell arrays
"""
struct CellArrayFromInner{A,B,T,N} <: CellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::CellArrayFromInner) = self.a

rightcellarray(self::CellArrayFromInner) = self.b

function computevals!(::CellArrayFromInner, a, b, v)
  v .= inner.(a,b)
end

"""
Lazy result of cellsum
"""
struct CellArrayFromCellSum{A,N,C,T} <: CellArrayFromUnaryOp{C,T,N}
  a::C
end

inputcellarray(self::CellArrayFromCellSum) = self.a

function computesize(self::CellArrayFromCellSum{A},asize) where A
  cellsumsize(asize,Val(A))
end

function computevals!(::CellArrayFromCellSum{A}, a, v) where A
  cellsumvals!(a,v,Val(A))
end

"""
Lazy result of cellnewaxis
"""
struct CellArrayFromCellNewAxis{A,C,T,N} <: CellArrayFromUnaryOp{C,T,N}
  a::C
end

inputcellarray(self::CellArrayFromCellNewAxis) = self.a

@generated function computesize(self::CellArrayFromCellNewAxis{A},asize::NTuple{M,Int}) where {A,M}
  @assert A <= M+1
  str = ["asize[$i]," for i in 1:M]
  insert!(str,A,"1,")
  Meta.parse("($(join(str)))")
end

@generated function computevals!(::CellArrayFromCellNewAxis{D}, A::AbstractArray{T,M}, B::AbstractArray{T,N}) where {D,T,M,N}
  @assert N == M + 1
  @assert D <= N
  quote
    @nloops $M a A begin
      @nexprs $N j->(b_j = j == $D ? 1 : a_{ j < $D ? j : j-1 } )
      (@nref $N B b) = @nref $M A a 
    end
  end    
end

