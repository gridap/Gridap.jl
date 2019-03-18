
export Point, MPoint
using StaticArrays

"""
Type representing a point of D dimensions
"""
const Point{D} = SVector{D,Float64} where D

"""
The mutable version of Point{D}
"""
const MPoint{D} = MVector{D,Float64} where D

