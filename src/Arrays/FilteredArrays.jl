
# This seems to be unused.
"""
"""
struct FilterMap <: Map end

function return_cache(k::FilterMap,f,a)
  # vals = testitem(a)
  vals = a
  T = eltype(eltype(a))
  r = zeros(T,length(vals))
  c = CachedArray(r)
end

function evaluate!(cache,k::FilterMap,f,a)
  c = cache
  vals = a
  filters = f
  @check size(vals) == size(filters) "Local arrays mismatch"
  setsize!(c,(sum(filters),))
  r = c.array
  i = 0
  for (val,filter) in zip(vals,filters)
    if filter
      i += 1
      r[i] = val
    end
  end
  r
end
