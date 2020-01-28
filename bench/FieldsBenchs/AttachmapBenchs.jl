module AttachmapBenchs

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: OtherMockBasis, MockBasis
using FillArrays
using Gridap.TensorValues

@inline function loop(a,cache)
  for i in eachindex(a)
    ai = getindex!(cache,a,i)
  end
end

function bench1(n)
  p1 = Point(2,2)
  p2 = Point(4,2)
  p3 = Point(1,3)
  p4 = Point(5,2)
  x = [p1,p2,p3,p4]
  np = length(x)
  v = 2.0
  d = 2
  ndof = 8
  wi = 3.0
  w = fill(wi,ndof)
  f = OtherMockBasis{d}(ndof)
  r = MockBasis{d}(v,ndof)
  l = n
  ax = fill(x,l)
  af = Fill(f,l)
  ar = Fill(r,l)
  aw = fill(w,l)
  aϕ = lincomb(af,aw)
  ab = attachmap(ar,aϕ)
  abx = evaluate(ab,ax)
  cabx = array_cache(abx)
  @time loop(abx,cabx)
  a∇bx = evaluate(∇(ab),ax)
  ca∇bx = array_cache(a∇bx)
  @time loop(a∇bx,ca∇bx)
end

for n in (1,1,10,1000,100000)
  @eval begin
    println("+++ runing suite for n = $($n) +++")
    bench1($n)
  end
end

end # module
