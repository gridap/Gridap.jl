module FieldValues

using TensorValues

export FieldValue
export Point
export VectorValue
export TensorValue
export MultiValue
export inner
export outer
export meas

"""
Type representing all possible field value types
"""
const FieldValue = Number

"""
Type representing a point of D dimensions with coordinates of type T
"""
const Point{D,T} = VectorValue{D,T}

end # module
