module IntegrateBenchs

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: IntKernel
using Gridap.Fields: OtherMockBasis, MockBasis, MockField
using FillArrays
using Gridap.TensorValues

@noinline function loop(a,cache)
  for i in eachindex(a)
    ai = getindex!(cache,a,i)
  end
  nothing
end

@inline function repeat(n,f,args...)
  for i in 1:n
    f(args...)
  end
  nothing
end

function bench1(n)

  fi = 2.0
  wi = 3.0
  ji = TensorValue(1.0,0.0,0.0,2.0)
  np = 4
  f = fill(fi,np)
  w = fill(wi,np)
  j = fill(ji,np)
  k = IntKernel()
  ck = kernel_cache(k,f,w,j)
  @time repeat(n,apply_kernel!,ck,k,f,w,j)

end

function bench2(n)

  fi = 2.0
  wi = 3.0
  ji = TensorValue(1.0,0.0,0.0,2.0)
  np = 4
  f = fill(fi,np,10)
  w = fill(wi,np)
  j = fill(ji,np)
  k = IntKernel()
  ck = kernel_cache(k,f,w,j)
  @time repeat(n,apply_kernel!,ck,k,f,w,j)

end

function bench3(n)

  np = 4
  p = Point(2,2)
  x = fill(p,np)

  ndof = 8
  d = 2
  ri = 5.0
  r = MockBasis{d}(ri,ndof)
  
  v = 3.0
  ndof = 8
  w = fill(1/np,np)
  c = fill(1.5,ndof)
  f = OtherMockBasis{d}(ndof)
  
  l = n
  ac = fill(c,l)
  af = Fill(f,l)
  aϕ = lincomb(af,ac)
  aj = ∇(aϕ)
  aw = Fill(w,l)
  ax = Fill(x,l)
  ar = Fill(r,l)
  ab = attachmap(ar,aϕ)

  fmass = field_array_operation(inner,ab,ab)
  mmass = integrate(fmass,ax,aw,aj)

  cmmass = array_cache(mmass)
  @time loop(mmass,cmmass)

  fstif = field_array_operation(inner,∇(ab),∇(ab))
  mstif = integrate(fstif,ax,aw,aj)

  cmstif = array_cache(mstif)
  @time loop(mstif,cmstif)

end

function bench4(n)

  fi = 2.0
  wi = 3.0
  ji = TensorValue(1.0,0.0,0.0,2.0)
  np = 4
  f = fill(fi,np)
  w = fill(wi,np)
  j = fill(ji,np)
  l = n
  af = fill(f,l)
  aw = fill(w,l)
  aj = fill(j,l)

  k = IntKernel()
  s = apply(k,af,aw,aj)

  cs = array_cache(s)
  @time loop(s,cs)

end

function bench5(n)

  fi = 2.0
  wi = 3.0
  ji = TensorValue(1.0,0.0,0.0,2.0)
  np = 4
  ndof = 10
  f = fill(fi,np,ndof)
  w = fill(wi,np)
  j = fill(ji,np)
  l = n
  af = fill(f,l)
  aw = fill(w,l)
  aj = fill(j,l)

  k = IntKernel()
  s = apply(k,af,aw,aj)

  cs = array_cache(s)
  @time loop(s,cs)
  @time repeat(n,getindex!,cs,s,1)

end

for n in (1,1,10,1000,100000)
  @eval begin
    println("+++ runing suite for n = $($n) +++")
    bench1($n)
    bench2($n)
    bench3($n)
    bench4($n)
    bench5($n)
  end
end

end #module
