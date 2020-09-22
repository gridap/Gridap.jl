abstract type Mapping <: GridapType end

return_cache(f,x...) = nothing

evaluate!(cache,f,x...) = @abstractmethod

# @fverdugo : TODO unify inference mechanism for Function and Mapping
return_type(f,x...) = typeof(testitem(f,x...))

function evaluate(f,x...)
  c = return_cache(f,x...)
  y = evaluate!(c,f,x...)
  y
end

(m::Mapping)(x...) = evaluate(m,x...)

function testitem(k,x...)
  cache = return_cache(k,x...)
  testitem!(cache,k,x...)
end

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

function return_caches(fs::Tuple,x...)
  _mapping_caches(x,fs...)
end

function _mapping_caches(x::Tuple,a,b...)
  ca = return_cache(a,x...)
  cb = return_caches(b,x...)
  (ca,cb...)
end

function _mapping_caches(x::Tuple,a)
  ca = return_cache(a,x...)
  (ca,)
end

@inline function evaluate!(cfs::Tuple,f::Tuple,x...)
  _evaluate_mappings!(cfs,x,f...)
end

function evaluate(fs::Tuple,x...)
cs = return_caches(fs,x...)
y = evaluate!(cs,fs,x...)
y
end

@inline function _evaluate_mappings!(cfs,x,f1,f...)
  cf1, cf = _split(cfs...)
  f1x = evaluate!(cf1,f1,x...)
  fx = evaluate!(cf,f,x...)
  (f1x,fx...)
end

@inline function _evaluate_mappings!(cfs,x,f1)
  cf1, = cfs
  f1x = evaluate!(cf1,f1,x...)
  (f1x,)
end

# @fverdugo : TODO this function is quite usefull. Better name and export it
@inline function _split(a,b...)
  (a,b)
end

function return_types(f::Tuple,x...)
  _mapping_return_types(x,f...)
end

function _mapping_return_types(x::Tuple,a,b...)
  Ta = return_type(a,x...)
  Tb = return_types(b,x...)
  (Ta,Tb...)
end

function _mapping_return_types(x::Tuple,a)
  Ta = return_type(a,x...)
  (Ta,)
end

function testitems(a,b...)
  va = testitem(a)
  vb = testitems(b...)
  (va,vb...)
end

function testitems(a)
  va = testitem(a)
  (va,)
end

# Extended Array interface

# @fverdugo: TODO better handling of Cartesian indices
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
  Ts = return_types(c.l,x...)
  return_type(c.k, testvalues(Ts...)...)
end

function return_cache(c::OperationMapping,x...)
  cl = return_caches(c.l,x...)
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
