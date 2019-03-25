export VectorValue, MVectorValue
export Point, MPoint
export FieldValue
export ScalarValue, VectorValue, TensorValue
export MVectorValue, MTensorValue
using StaticArrays

"""
Type representing a scalar value
"""
const ScalarValue = Float64

"""
Type representing a vector value of dimension `D`
"""
const VectorValue{D} = SVector{D,Float64} where D

"""
Type representing a tensor value of dimension `D`
"""
const TensorValue{D,DD} = SMatrix{D,D,Float64,DD} where D
# @santiagobadia : Any way to eliminate DD?

"""
Mutable version of `VectorValue{D}`
"""
const MVectorValue{D} = MVector{D,Float64} where D

"""
Mutable version of `TensorValue{D,DD}`
"""
const MTensorValue{D,DD} = MMatrix{D,D,Float64,DD} where D
# @santiagobadia : Any way to eliminate DD?

"""
Type representing all possible field value types
"""
const FieldValue = Union{ScalarValue, VectorValue, TensorValue, MVectorValue, MTensorValue}
# @santiagobadia : Should we make a difference between SFieldValue and MFieldValue?

"""
Type representing a point of D dimensions
"""
const Point{D} = SVector{D,Float64} where D

"""
The mutable version of Point{D}
"""
const MPoint{D} = MVector{D,Float64} where D

gradient(::Type{Float64},::Val{D}) where D = VectorValue{D}

gradient(::Type{VectorValue{D}},::Val{Z}) where {D,Z} = SMatrix{Z,D,Float64,Z*D}

outer(a::T,b::T) where T <: Number = a*b

inner(a::T,b::T) where T <: Number = a*b

outer(a::T,b::SVector{D,T}) where {T <: Number,D} = a*b

function outer(a::SVector{D,T},b::SVector{Z,T}) where {D,Z,T}
  if D==2 && Z==2
    SMatrix{2,2,T,4}( a[1]*b[1], a[2]*b[1], a[1]*b[2], a[2]*b[2] )
  else
    @notimplemented
  end
end

function inner(a::SVector{D,T},b::SVector{D,T}) where {D,T}
  v = 0.0
  @inbounds for i=1:D
    v += a[i]*b[i]
  end
  v
end
