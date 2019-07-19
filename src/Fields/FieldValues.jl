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
export normalvec

import TensorValues: meas


"""
Type representing all possible field value types
"""
const FieldValue = Number

"""
Type representing a point of D dimensions with coordinates of type T
"""
const Point{D,T} = VectorValue{D,T}

function normalvec(v::MultiValue{Tuple{1,2}})
  n1 = v[1,2]
  n2 = -1*v[1,1]
  VectorValue(n1,n2)
end

function normalvec(v::MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  VectorValue(n1,n2,n3)
end

function meas(v::MultiValue{Tuple{1,2}})
  n = normalvec(v)
  sqrt(n*n)
end

function meas(v::MultiValue{Tuple{2,3}})
  n = normalvec(v)
  sqrt(n*n)
end

end # module
