function fill_entries!(J::AbstractArray,v)
  J .= convert(eltype(J),v)
end
