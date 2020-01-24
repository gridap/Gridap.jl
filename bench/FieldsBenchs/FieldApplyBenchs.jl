module FieldApplyBenchs

using Gridap.Arrays
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
  g = apply_kernel_to_field(bcast(+),f,v)
  cg = field_cache(g,x)
  @time repeat(n,evaluate!,cg,g,x)
  ∇g = gradient(g)
  ∇cg = field_cache(∇g,x)
  @time repeat(n,evaluate!,∇cg,∇g,x)
  h = apply_kernel_to_field(bcast(+),f,f)
  ch = field_cache(h,x)
  @time repeat(n,evaluate!,ch,h,x)
end

#using Gridap.Inference
#
#function bench3(n)
#  np = 4
#  p = Point(1,2)
#  x = fill(p,np)
#  v = 3.0
#  d = 2
#  f = MockField{d}(v)
#  k = bcast(+)
#  @time repeat(n,apply_kernel_to_field,k,f,v)
#end
#
#using Gridap.Fields: Valued
#
#function bench4(n)
#  np = 4
#  p = Point(1,2)
#  x = fill(p,np)
#  v = 3.0
#  ndof = 8
#  w = fill(v,ndof)
#  d = 2
#  f = MockBasis{d}(v,ndof)
#  k = bcast(+)
#  fx = evaluate(f,x)
#  T = kernel_return_type(k,fx,fx)
#  @time repeat(n,apply_kernel_to_field,T,k,f,f)
#  @time repeat(n,apply_kernel_to_field,T,k,f,w)
#  @time repeat(n,apply_kernel_to_field,T,k,w,f)
#  vk = Valued(k,f,w)
#  cvk = kernel_cache(vk,f,w)
#  @time repeat(n,apply_kernel!,cvk,vk,f,w)
#end

for n in (1,1,10,1000,100000)
  @eval begin
    println("+++ runing suite for n = $($n) +++")
    bench1($n)
    #bench3($n)
    #bench4($n)
  end
end

end # module
