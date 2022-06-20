struct KeyToVal{T} <: Map
  key_to_val::T
end


function return_cache(m::KeyToVal,keys)
  K = eltype(keys)
  length(keys) > 0 ? _key = keys[1] : _key = testvalue(K)
  V = typeof(m.key_to_val(_key))
  d = Dict{K,V}()
end

function evaluate!(cache,m::KeyToVal,key)
  if haskey(cache,key)
    return get(cache,key,nothing)
  else
    val = m.key_to_val(key)
    cache[key] = val
    return val
  end
end
