function fill_entries!(J::AbstractArray,v::Number)
  J .= convert(eltype(J),v)
end
