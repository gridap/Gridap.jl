module FieldArraysBenchs

using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: MockField, MockBasis
using FillArrays
using Gridap.TensorValues

@noinline function loop(a,cache)
  for i in eachindex(a)
    ai = getindex!(cache,a,i)
  end
end

@noinline function loop_and_evaluate(ca,cai,cx,a,x)
  for i in eachindex(a)
    ai = getindex!(ca,a,i)
    xi = getindex!(cx,x,i)
    vi = evaluate!(cai,ai,xi)
  end
end

function bench1(n)
  np = 4
  p = Point(1,2)
  x = fill(p,np)
  v = 3.0
  d = 2
  l = n
  f = MockField{d}(v)
  af = fill(f,l)
  ax = fill(x,l)
  afx = evaluate(af,ax)
  cafx = array_cache(afx)
  @time loop(afx,cafx)
  ag = apply_to_field_array(bcast(+),af,af)
  cag = array_cache(ag)
  @time loop(ag,cag)
  agx = evaluate(ag,ax)
  cagx = array_cache(agx)
  @time loop(agx,cagx)
  g = testitem(ag)
  cg = field_cache(g,x)
  cax = array_cache(ax)
  @time loop_and_evaluate(cag,cg,cax,ag,ax)
end

function bench2(n)
  np = 4
  p = Point(1,2)
  x = fill(p,np)
  v = 3.0
  d = 2
  l = n
  f = MockField{d}(v)
  af = Fill(f,l)
  ax = fill(x,l)
  afx = evaluate(af,ax)
  cafx = array_cache(afx)
  @time loop(afx,cafx)
  _ag = apply_to_field_array(bcast(+),af,af)
  ag = gradient(_ag)
  cag = array_cache(ag)
  @time loop(ag,cag)
  agx = evaluate(ag,ax)
  cagx = array_cache(agx)
  @time loop(agx,cagx)
  g = testitem(ag)
  cg = field_cache(g,x)
  cax = array_cache(ax)
  @time loop_and_evaluate(cag,cg,cax,ag,ax)
end

function bench2a(n)
  np = 4
  p = Point(1,2)
  x = fill(p,np)
  v = 3.0
  d = 2
  l = n
  f = MockField{d}(v)
  af = Fill(f,l)
  av = fill(v,l)
  ax = fill(x,l)
  afx = evaluate(af,ax)
  cafx = array_cache(afx)
  @time loop(afx,cafx)
  _ag = apply_to_field_array(bcast(+),af,av)
  ag = gradient(_ag)
  cag = array_cache(ag)
  @time loop(ag,cag)
  agx = evaluate(ag,ax)
  cagx = array_cache(agx)
  @time loop(agx,cagx)
  g = testitem(ag)
  cg = field_cache(g,x)
  cax = array_cache(ax)
  @time loop_and_evaluate(cag,cg,cax,ag,ax)
end

function bench3(n)
  np = 4
  p = Point(1,2)
  x = fill(p,np)
  v = 3.0
  d = 2
  ndof = 8
  wi = 3.0
  w = fill(wi,ndof)
  l = n
  f = MockBasis{d}(v,ndof)
  af = Fill(f,l)
  ax = fill(x,l)
  aw = Fill(w,l)
  afx = evaluate(af,ax)
  cafx = array_cache(afx)
  @time loop(afx,cafx)
  ag = apply_to_field_array(bcast(+),af,aw)
  cag = array_cache(ag)
  @time loop(ag,cag)
  agx = evaluate(ag,ax)
  cagx = array_cache(agx)
  @time loop(agx,cagx)
  g = testitem(ag)
  cg = field_cache(g,x)
  cax = array_cache(ax)
  @time loop_and_evaluate(cag,cg,cax,ag,ax)
end

for n in (1,1,10,1000,100000)
  @eval begin
    println("+++ runing suite for n = $($n) +++")
    bench1($n)
    bench2($n)
    bench2a($n)
    bench3($n)
  end
end

end # module
