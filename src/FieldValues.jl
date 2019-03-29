module FieldValues

export FieldValue
export Point, MPoint
export ScalarValue, VectorValue, TensorValue
export MVectorValue, MTensorValue
export inner, outer

using StaticArrays: SVector, MVector, SMatrix, MMatrix

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
const TensorValue{D,DD} = SMatrix{D,D,Float64,DD} where {D,DD}

"""
Mutable version of `VectorValue{D}`
"""
const MVectorValue{D} = MVector{D,Float64} where D

"""
Mutable version of `TensorValue{D,DD}`
"""
const MTensorValue{D,DD} = MMatrix{D,D,Float64,DD} where {D,DD}

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

# Constructors

@generated function VectorValue(args::Vararg{Float64,D}) where D
  :( VectorValue{$D}(args...)  )
end

@generated function MVectorValue(args::Vararg{Float64,D}) where D
  :( MVectorValue{$D}(args...)  )
end

@generated function TensorValue(args::Vararg{Float64,DD}) where DD
  SQ = sqrt(DD)
  D = ceil(Int,SQ)
  @assert D == SQ
  :( TensorValue{$D,$DD}(args...)  )
end

@generated function MTensorValue(args::Vararg{Float64,DD}) where DD
  SQ = sqrt(DD)
  D = ceil(Int,SQ)
  @assert D == SQ
  :( MTensorValue{$D,$DD}(args...)  )
end

Point(x) = VectorValue(x)

MPoint(x) = MVectorValue(x)

# Operations

outer(a::T,b::T) where T <: Number = a*b

outer(::Type{T},::Type{T}) where T <: Number = T

outer(a::T,b::SVector{D,T}) where {T <: Number,D} = a*b

outer(::Type{T},::Type{SVector{D,T}}) where {T <: Number,D} = SVector{D,T}

outer(b::SVector{D,T},a::T) where {T <: Number,D} = a*b

outer(::Type{SVector{D,T}},::Type{T}) where {T <: Number,D} = SVector{D,T}

outer(a::T,b::SMatrix{D,E,T,DE}) where {T <: Number,D,E,DE} = a*b

outer(::Type{T},::Type{SMatrix{D,E,T,DE}}) where {T <: Number,D,E,DE} = SMatrix{D,E,T,DE}

outer(b::SMatrix{D,E,T,DE},a::T) where {T <: Number,D,E,DE} = a*b

outer(::Type{SMatrix{D,E,T,DE}},::Type{T}) where {T <: Number,D,E,DE} = SMatrix{D,E,T,DE}

@generated function outer(a::SVector{D,T},b::SVector{Z,T}) where {D,Z,T}
  str = join(["a[$i]*b[$j], " for j in 1:Z for i in 1:D])
  Meta.parse("SMatrix{$D,$Z,Float64,$(D*Z)}($str)")
end

@generated function outer(::Type{SVector{D,T}},::Type{SVector{Z,T}}) where {D,Z,T}
  Meta.parse("SMatrix{$D,$Z,$T,$(D*Z)}")
end

inner(a::T,b::T) where T <: Number = a*b

inner(::Type{T},::Type{T}) where T <: Number = T

@generated function inner(a::SVector{D,T},b::SVector{D,T}) where {D,T}
  str = join([" a[$i]*b[$i] +" for i in 1:D])
  Meta.parse(str[1:(end-1)])
end

inner(::Type{SVector{D,T}},::Type{SVector{D,T}}) where {D,T} = T

@generated function inner(a::SMatrix{D,Z,T,DZ},b::SMatrix{D,Z,T,DZ}) where {D,Z,T,DZ}
  str = join([" a[$i,$j]*b[$i,$j] +" for j in 1:Z for i in 1:D])
  Meta.parse(str[1:(end-1)])
end

inner(::Type{SMatrix{D,Z,T,DZ}},::Type{SMatrix{D,Z,T,DZ}}) where {D,Z,T,DZ} = T

end # module FieldValues
