module CachedValues

export CachedValue

mutable struct CachedValue{F}
  value::F
end

end # module
