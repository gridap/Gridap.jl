"""
    struct KeyToValMap{T<:Function} <: Map

Map for lazily filling a `Dict` the outputs of the function `T` over an array of inputs.
"""
struct KeyToValMap{T<:Function} <: Map
  key_to_val::T
end

function return_cache(m::KeyToValMap,key)
  K = typeof(key)
  V = typeof(m.key_to_val(key))
  d = Dict{K,V}()
end

function evaluate!(cache,m::KeyToValMap,key)
  if haskey(cache,key)
    return get(cache,key,nothing)
  else
    val = m.key_to_val(key)
    cache[key] = val
    return val
  end
end
