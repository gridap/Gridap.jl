module ArraysBenchs

using Gridap.Arrays
using FillArrays

@inline function loop(a,cache)
  for i in eachindex(a)
    ai = getindex!(cache,a,i)
  end
end

@inline function loop_and_apply(ca,cai,cx,a,x...)
  for i in eachindex(a)
    ai = getindex!(ca,a,i)
    xi = getitems!(cx,x,i...)
    vi = apply_kernel!(cai,ai,xi...)
  end
end

function bench1(n)
  a = rand(n)
  k = -
  c = apply(k,a)
  cache = array_cache(c)
  @time loop(c,cache)
end

function bench2(n)
  a = rand(n)
  b = rand(n)
  k = -
  c = apply(k,a,b)
  cache = array_cache(c)
  @time loop(c,cache)
end

function bench3(n)
  a = fill(rand(2,3),n)
  b = rand(n)
  c = apply(bcast(-),a,b)
  cache = array_cache(c)
  @time loop(c,cache)
end

function bench4(n)
  a = fill(rand(2,3),n)
  b = rand(n)
  c = apply(bcast(-),a,b)
  d = apply(bcast(+),a,c)
  e = apply(bcast(*),d,c)
  cache = array_cache(e)
  @time loop(e,cache)
end

function bench7(n)
  k = +
  a = fill(k,n)
  x = rand(n)
  y = rand(n)
  v = apply(a,x,y)
  cache = array_cache(v)
  @time loop(v,cache)
end

function bench8(n)
  a = fill(bcast(+),n)
  x = [rand(2,3) for i in 1:n]
  y = [rand(1,3) for i in 1:n]
  v = apply(a,x,y)
  cache = array_cache(v)
  @time loop(v,cache)
end

function bench9(n)
  a = Fill(bcast(+),n)
  x = [rand(mod(i-1,3)+1,3) for i in 1:n]
  y = [rand(1,3) for i in 1:n]
  v = apply(a,x,y)
  cache = array_cache(v)
  @time loop(v,cache)
end

for n in (1,1,10,1000,100000)
  @eval begin
    println("+++ running suite for n = $($n) +++")
    bench1($n)
    bench2($n)
    bench3($n)
    bench4($n)
    bench7($n)
    bench8($n)
    bench9($n)
  end
end

end # module
