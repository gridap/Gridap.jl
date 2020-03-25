
function reinterpret(a::Array{VectorValue{D,T}}) where {D,T}
  b = reinterpret(T,a)
  sa = size(a)
  t = size(VectorValue{D,T})
  s = (t...,sa...)
  reshape(b,s)
end

function reinterpret(a::Array{TensorValue{D1,D2,T,L}}) where {D1,D2,T,L}
  b = reinterpret(T,a)
  sa = size(a)
  t = size(TensorValue{D1,D2,T,L})
  s = (t...,sa...)
  reshape(b,s)
end

