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
export trace
export tr
export symmetic_part

import TensorValues: meas
import Base: adjoint


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

@generated function trace(v::TensorValue{D}) where D
  str = join([" v.array.data[$i+$((i-1)*D)] +" for i in 1:D ])
  Meta.parse(str[1:(end-1)])
end

const tr = trace

@generated function symmetic_part(v::TensorValue{D}) where D
  str = "("
  for j in 1:D
    for i in 1:D
      str *= "0.5*v.array.data[$i+$((j-1)*D)] + 0.5*v.array.data[$j+$((i-1)*D)], "
    end
  end
  str *= ")"
  Meta.parse("TensorValue($str)")
end

function adjoint(v::TensorValue)
  t = adjoint(v.array)
  TensorValue(t)
end

end # module
