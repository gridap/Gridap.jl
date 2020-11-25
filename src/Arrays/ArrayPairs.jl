"""
"""
function pair_arrays(a::AbstractArray,b::AbstractArray)
  lazy_map(tuple,a,b)
end

"""
"""
function unpair_arrays(pair::AbstractArray{<:Tuple{<:Any,<:Any}})
  a = lazy_map(first,pair)
  b = lazy_map(second,pair)
  a, b
end

second(x) = x[2]

function unpair_arrays(pair::LazyArray{<:Fill{typeof(tuple)}})
  pair.f
end

