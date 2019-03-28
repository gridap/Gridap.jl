
function Base.:+(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  CellArrayFromSum{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:-(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  CellArrayFromSub{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:*(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  CellArrayFromMul{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:/(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  CellArrayFromDiv{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:(==)(a::CellArray{T,N},b::CellArray{T,N}) where {T,N}
  length(a) != length(b) && return false
  cellsize(a) != cellsize(b) && return false
  if N != 1; @notimplemented end
  for ((ai,ais),(bi,bis)) in zip(a,b)
    ais != bis && return false
    for j in 1:ais[1]
      ai[j] != bi[j] && return false
    end
  end
  return true
end

"""
Assumes that det is defined for instances of T
and that the result is Float64
"""
function LinearAlgebra.det(self::CellArray{T,N}) where {T,N}
  CellArrayFromDet{typeof(self),Float64,N}(self)
end

"""
Assumes that inv is defined for instances of T
"""
function LinearAlgebra.inv(self::CellArray{T,N}) where {T,N}
  CellArrayFromInv{typeof(self),T,N}(self)
end

"""
Type that stores the lazy result of evaluating the determinant
of each element in a CellArray
"""
struct CellArrayFromDet{C,T,N} <: CellArrayFromElemUnaryOp{C,T,N}
  a::C
end

inputcellarray(self::CellArrayFromDet) = self.a

function computevals!(::CellArrayFromDet, a, asize, v, vsize)
  v .= det.(a)
end

"""
Type that stores the lazy result of evaluating the inverse of
of each element in a CellArray
"""
struct CellArrayFromInv{C,T,N} <: CellArrayFromElemUnaryOp{C,T,N}
  a::C
end

inputcellarray(self::CellArrayFromInv) = self.a

function computevals!(::CellArrayFromInv, a, asize, v, vsize)
  v .= inv.(a)
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

function computevals!(::CellArrayFromSum, a, asize, b, bsize, v, vsize)
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

function computevals!(::CellArrayFromSub, a, asize, b, bsize, v, vsize)
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

function computevals!(::CellArrayFromMul, a, asize, b, bsize, v, vsize)
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

function computevals!(::CellArrayFromDiv, a, asize, b, bsize, v, vsize)
  v .= a ./ b
end
