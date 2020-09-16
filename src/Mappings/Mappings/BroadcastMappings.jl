struct BroadcastMapping{F<:Function} <: Mapping
  f::F
end

function return_type(f::BroadcastMapping,x::Number...)
  Ts = map(typeof,x)
  return_type(f.f,Ts...)
end

function return_cache(f::BroadcastMapping,x::Number...)
  nothing
end

@inline function evaluate!(::Nothing,f::BroadcastMapping,x::Number...)
  f.f(x...)
end

function return_type(f::BroadcastMapping,x::NumberOrArray...)
  typeof(return_cache(f,x...).array)
end

function return_cache(f::BroadcastMapping,x::NumberOrArray...)
  s = _sizes(x...)
  bs = Base.Broadcast.broadcast_shape(s...)
  Te = map(numbertype,x)
  T = return_type(f.f,Te...)
  N = length(bs)
  r = Array{T,N}(undef,bs)
  ri = testvalue(T)
  fill!(r,ri)
  cache = CachedArray(r)
   _prepare_cache(cache,x...)
end

numbertype(a::AbstractArray) = eltype(a)

numbertype(a::Number) = typeof(a)

@inline function evaluate!(cache,f::BroadcastMapping,x::NumberOrArray...)
  r = _prepare_cache(cache,x...)
  a = r.array
  broadcast!(f.f,a,x...)
  a
end

@inline function _prepare_cache(c,x...)
  s = _sizes(x...)
  bs = Base.Broadcast.broadcast_shape(s...)
  if bs != size(c)
    setsize!(c,bs)
  end
  c
end

# TODO use map
@inline function _sizes(a,x...)
  (_size(a), _sizes(x...)...)
end

@inline function _sizes(a)
  (_size(a),)
end

@inline _size(a) = size(a)
@inline _size(a::Number) = (1,)

function _checks(a,b)
  @assert size(a) == size(b) "Sizes must agree."
  nothing
end
