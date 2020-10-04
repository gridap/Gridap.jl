abstract type Mapping <: GridapType end

return_cache(f,x...) = nothing

evaluate!(cache,f,x...) = @abstractmethod

return_type(f,x...) = typeof(testitem(f,x...))

function evaluate(f,x...)
  c = return_cache(f,x...)
  y = evaluate!(c,f,x...)
  y
end

(m::Mapping)(x...) = evaluate(m,x...)

@inline testitem(k,x...) = evaluate(k,x...)

@inline function testitem!(cache,k,x...)
  evaluate!(cache,k,x...)
end

# Default implementation for Function

evaluate!(cache,f::Function,x...) = f(x...)

# Number or Array implementation

const NumberOrArray = Union{Number,AbstractArray{<:Number}}

evaluate!(cache,f::NumberOrArray,x...) = f
return_type(f::NumberOrArray,x...) = typeof(f)

# Testing the interface

function test_mapping(f,x::Tuple,y,cmp=(==))
  z = evaluate(f,x...)
  @test cmp(z,y)
  @test typeof(z) == return_type(f,x...)
  cache = return_cache(f,x...)
  z = evaluate!(cache,f,x...)
  @test cmp(z,y)
  z = evaluate!(cache,f,x...)
  @test cmp(z,y)
  z = testitem!(cache,f,x...)
  @test cmp(typeof(z),typeof(y))
end

# Work with several Mapping objects

@inline function evaluate!(cfs::Tuple,f::Tuple,x...)
  map((c,fi) -> evaluate!(c,fi,x...),cfs,f)
  # _evaluate_mappings!(cfs,x,f...)
end

function evaluate(fs::Tuple,x...)
cs = map(fi -> return_cache(fi,x...),fs)
y = evaluate!(cs,fs,x...)
y
end

# Extended Array interface

function return_cache(a::AbstractArray)
  i = testitem(eachindex(a))
  return_cache(a,Tuple(i)...)
end

function return_cache(a::AbstractArray,i...)
  nothing
end

# Broadcast Functions

struct BroadcastMapping{F} <: Mapping
  f::F
end

@inline function evaluate!(cache,f::BroadcastMapping,x...)
  r = _prepare_cache(cache,x...)
  a = r.array
  broadcast!(f.f,a,x...)
  a
end

function evaluate!(cache,b::BroadcastMapping,args::Number...)
  b.f(args...)
end

function return_type(f::BroadcastMapping,x::Number...)
  Ts = map(typeof,x)
  return_type(f.f,Ts...)
end

function return_type(f::BroadcastMapping,x::AbstractArray...)
  typeof(return_cache(f,x...).array)
end

function return_cache(f::BroadcastMapping,x::Number...)
  nothing
end

function return_cache(f::BroadcastMapping,x...)
  s = _size.(x)
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

@inline function _prepare_cache(c,x...)
  s = _size.(x)
  bs = Base.Broadcast.broadcast_shape(s...)
  if bs != size(c)
    setsize!(c,bs)
  end
  c
end

@inline _size(a) = size(a)
@inline _size(a::Number) = (1,)

function _checks(a,b)
  @assert size(a) == size(b) "Sizes must agree."
  nothing
end

# OperationMappings

@inline function composition(k,l...)
  OperationMapping(k,l)
end

struct OperationMapping{K,L} <: Mapping
  k::K
  l::L
  @inline function OperationMapping(k,l)
    new{typeof(k),typeof(l)}(k,l)
  end
end

function return_type(c::OperationMapping,x...)
  Ts = map(fi -> return_type(fi,x...),c.l)
  # Ts = return_types(c.l,x...)
  return_type(c.k, testvalues(Ts...)...)
end

function return_cache(c::OperationMapping,x...)
  cl = map(fi -> return_cache(fi,x...),c.l)
  # cl = return_caches(c.l,x...)
  lx = evaluate!(cl,c.l,x...)
  ck = return_cache(c.k,lx...)
  (ck,cl)
end

@inline function evaluate!(cache,c::OperationMapping,x...)
  ck, cf = cache
  lx = evaluate!(cf,c.l,x...)
  evaluate!(ck,c.k,lx...)
end

# Operations

struct Operation{T} <: Mapping
  op::T
end

evaluate!(cache,op::Operation,x...) = OperationMapping(op.op,x)

(op::Operation)(x...) = evaluate!(nothing,op,x...)

operation(x) = Operation(x)
