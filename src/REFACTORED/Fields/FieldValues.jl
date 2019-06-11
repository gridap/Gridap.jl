module FieldValues

using TensorValues

export FieldValue
export Point

"""
Type representing all possible field value types
"""
const FieldValue = Number

"""
Type representing a point of D dimensions with coordinates of type T
"""
const Point{D,T} = VectorValue{D,T}

end # module
