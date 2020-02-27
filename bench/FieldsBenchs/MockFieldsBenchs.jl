module MockFieldsBenchs

using Gridap.Fields
using Gridap.Fields: MockField, MockBasis

@inline function repeat(n,f,args...)
  for i in 1:n
    f(args...)
  end
  nothing
end

function bench1(n)
  np = 4
  p = Point(1,2)
  x = fill(p,np)
  v = 3.0
  d = 2
  f = MockField{d}(v)
  cache = field_cache(f,x)
  @time repeat(n,evaluate!,cache,f,x)
  ∇f = gradient(f)
  ∇cache = field_cache(∇f,x)
  @time repeat(n,evaluate!,∇cache,∇f,x)
end

function bench2(n)
  np = 4
  p = Point(1,2)
  x = fill(p,np)
  v = 3.0
  d = 2
  ndof = 8
  f = MockBasis{d}(v,ndof)
  cache = field_cache(f,x)
  @time repeat(n,evaluate!,cache,f,x)
  ∇f = gradient(f)
  ∇cache = field_cache(∇f,x)
  @time repeat(n,evaluate!,∇cache,∇f,x)
end

for n in (1,1,10,1000,100000)
  @eval begin
    println("+++ running suite for n = $($n) +++")
    bench1($n)
    bench2($n)
  end
end

end # module
