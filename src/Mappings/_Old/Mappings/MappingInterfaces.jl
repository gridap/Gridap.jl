abstract type Mapping <: GridapType end

return_type(f,x...) = typeof(testitem(f,x...))

return_cache(f,x...) = nothing

evaluate!(cache,f,x...) = @abstractmethod

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

# Function implementation

evaluate!(cache,f::Function,x...) = f(x...)

# @fverdugo : TODO This ad-hoc definition should be not needded
return_type(f::Function,x...) = return_type(f,map(typeof,x))

return_cache(f::Function,x...) = nothing

evaluate(f::Function,x...) = f(x...)

evaluate(T::Type,f::Function,x...) = f(x...)

testitem(f::Function,x...) = f(x...)

testitem!(cache,f::Function,x...) = f(x...)

# Number or Array implementation

const NumberOrArray = Union{Number,AbstractArray{<:Number}}

evaluate!(cache,f::NumberOrArray,x...) = f

return_cache(f::NumberOrArray,x...) = nothing

return_type(f::NumberOrArray,x...) = typeof(f)

evaluate(f::NumberOrArray,x...) = f

evaluate(T::Type,f::NumberOrArray,x...) = f

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
