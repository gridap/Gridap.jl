module MonomialBasesBenchs

using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials

@inline function repeat(n,f,args...)
  for i in 1:n
    f(args...)
  end
  nothing
end

function bench1(n)

  xi = Point(2,3)
  np = 5
  x = fill(xi,np)

  orders = (1,2)
  V = Float64
  b = MonomialBasis{2}(V,orders)

  cb = field_cache(b,x)
  @time repeat(n,evaluate_field!,cb,b,x)

  ∇b = ∇(b)
  c∇b = field_cache(∇b,x)
  @time repeat(n,evaluate_field!,c∇b,∇b,x)

end

function bench2(n)

  xi = Point(2,3)
  np = 5
  x = fill(xi,np)

  orders = (1,2)
  V = VectorValue{2,Float64}
  b = MonomialBasis{2}(V,orders)

  cb = field_cache(b,x)
  @time repeat(n,evaluate_field!,cb,b,x)

  ∇b = ∇(b)
  c∇b = field_cache(∇b,x)
  @time repeat(n,evaluate_field!,c∇b,∇b,x)

end

for n in (1,1,10,1000,100000)
  @eval begin
    println("+++ running suite for n = $($n) +++")
    bench1($n)
    bench2($n)
  end
end

end # module
