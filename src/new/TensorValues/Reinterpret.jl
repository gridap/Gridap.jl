
function reinterpret(a::Array{MultiValue{S,T,N,L}}) where {S,T,N,L}
  b = reinterpret(T,a)
  sa = size(a)
  sv = Size(S)
  t = _Size_to_tuple(sv)
  s = (t...,sa...)
  reshape(b,s)
end

_Size_to_tuple(::Size{t}) where t = t
