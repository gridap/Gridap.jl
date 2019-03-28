
function Base.:+(a::OtherCellArray{T,N},b::OtherCellArray{T,N}) where {T,N}
  OtherCellArrayFromSum{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:-(a::OtherCellArray{T,N},b::OtherCellArray{T,N}) where {T,N}
  OtherCellArrayFromSub{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:*(a::OtherCellArray{T,N},b::OtherCellArray{T,N}) where {T,N}
  OtherCellArrayFromMul{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:/(a::OtherCellArray{T,N},b::OtherCellArray{T,N}) where {T,N}
  OtherCellArrayFromDiv{typeof(a),typeof(b),T,N}(a,b)
end

function Base.:(==)(a::OtherCellArray{T,N},b::OtherCellArray{T,N}) where {T,N}
  length(a) != length(b) && return false
  maxsize(a) != maxsize(b) && return false
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
function LinearAlgebra.det(self::OtherCellArray{T,N}) where {T,N}
  OtherCellArrayFromDet{typeof(self),Float64,N}(self)
end

"""
Assumes that inv is defined for instances of T
"""
function LinearAlgebra.inv(self::OtherCellArray{T,N}) where {T,N}
  OtherCellArrayFromInv{typeof(self),T,N}(self)
end

"""
Type that stores the lazy result of evaluating the determinant
of each element in a CellArray
"""
struct OtherCellArrayFromDet{C,T,N} <: OtherCellArrayFromElemUnaryOp{C,T,N}
  a::C
end

inputcellarray(self::OtherCellArrayFromDet) = self.a

function computevals!(::OtherCellArrayFromDet, a, asize, v, vsize)
  v .= det.(a)
end

"""
Type that stores the lazy result of evaluating the inverse of
of each element in a CellArray
"""
struct OtherCellArrayFromInv{C,T,N} <: OtherCellArrayFromElemUnaryOp{C,T,N}
  a::C
end

inputcellarray(self::OtherCellArrayFromInv) = self.a

function computevals!(::OtherCellArrayFromInv, a, asize, v, vsize)
  v .= inv.(a)
end

"""
Lazy sum of two cell arrays
"""
struct OtherCellArrayFromSum{A,B,T,N} <: OtherCellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::OtherCellArrayFromSum) = self.a

rightcellarray(self::OtherCellArrayFromSum) = self.b

function computevals!(::OtherCellArrayFromSum, a, asize, b, bsize, v, vsize)
  v .= a .+ b
end

"""
Lazy subtraction of two cell arrays
"""
struct OtherCellArrayFromSub{A,B,T,N} <: OtherCellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::OtherCellArrayFromSub) = self.a

rightcellarray(self::OtherCellArrayFromSub) = self.b

function computevals!(::OtherCellArrayFromSub, a, asize, b, bsize, v, vsize)
  v .= a .- b
end

"""
Lazy multiplication of two cell arrays
"""
struct OtherCellArrayFromMul{A,B,T,N} <: OtherCellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::OtherCellArrayFromMul) = self.a

rightcellarray(self::OtherCellArrayFromMul) = self.b

function computevals!(::OtherCellArrayFromMul, a, asize, b, bsize, v, vsize)
  v .= a .* b
end

"""
Lazy division of two cell arrays
"""
struct OtherCellArrayFromDiv{A,B,T,N} <: OtherCellArrayFromElemBinaryOp{A,B,T,N}
  a::A
  b::B
end

leftcellarray(self::OtherCellArrayFromDiv) = self.a

rightcellarray(self::OtherCellArrayFromDiv) = self.b

function computevals!(::OtherCellArrayFromDiv, a, asize, b, bsize, v, vsize)
  v .= a ./ b
end
