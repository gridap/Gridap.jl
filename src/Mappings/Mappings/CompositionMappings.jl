@inline function composition(k,l...)
  CompositionMapping(k,l...)
end

struct CompositionMapping{K,L} <: Mapping
  k::K
  l::L
  @inline function CompositionMapping(k,l...)
    new{typeof(k),typeof(l)}(k,l)
  end
end

function return_type(c::CompositionMapping,x)
  Ts = return_types(c.l,x)
  return_type(c.k, testvalues(Ts...)...)
end

function return_cache(c::CompositionMapping,x)
  cl = return_caches(c.l,x)
  lx = evaluate!(cl,c.l,x)
  ck = return_cache(c.k,lx...)
  (ck,cl)
end

@inline function evaluate!(cache,c::CompositionMapping,x)
  ck, cf = cache
  lx = evaluate!(cf,c.l,x)
  evaluate!(ck,c.k,lx...)
end
