
struct FilterMap <: Map end

function return_cache(::FilterMap,f,a)
  T = eltype(eltype(a))
  CachedArray(zeros(T,length(a)))
end

function evaluate!(cache,::FilterMap,f,a)
  @check size(a) == size(f) "Local arrays mismatch"
  setsize!(cache,(sum(f),))
  r = cache.array
  i = 0
  for (fi,ai) in zip(f,a)
    if fi
      i += 1
      r[i] = ai
    end
  end
  r
end
