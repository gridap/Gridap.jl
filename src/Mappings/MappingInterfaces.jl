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

function testitem(k::Mapping,x...)
  cache = return_cache(k,x...)
  testitem!(cache,k,x...)
end

@inline function testitem!(cache,k::Mapping,x...)
  evaluate!(cache,k,x...)
end

# Default implementation for Function

evaluate!(cache,f::Function,x...) = f(x...)
# @fverdugo : TODO This ad-hoc definition should be not needded
return_type(f::Function,x...) = return_type(f,map(typeof,x))

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
