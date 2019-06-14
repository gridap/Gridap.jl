module CachedStructFields

export CachedStructField

mutable struct CachedStructField{F}
  value::F
end

end # module
