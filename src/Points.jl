
export VectorValue, MVectorValue
export Point, MPoint
using StaticArrays

"""
Type representing a vector value of D components
"""
const VectorValue{D} = SVector{D,Float64} where D

"""
Mutable version of VectorValue
"""
const MVectorValue{D} = MVector{D,Float64} where D

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

