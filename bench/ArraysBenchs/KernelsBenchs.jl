module KernelsBenchs

using Gridap.Arrays

@inline function repeat(n,f,args...)
  for i in 1:n
    f(args...)
  end
  nothing
end

function bench1(n)
  a = 1
  b = 2
  k = +
  cache = kernel_cache(k,a,b)
  @time repeat(n,apply_kernel!,cache,k,a,b)
end

function bench2(n)
  f = bcast(+)
  a = rand(3,2)
  b = 3
  cache = kernel_cache(f,a,b)
  @time repeat(n,apply_kernel!,cache,f,a,b)
end

function bench3(n)
  f = elem(+)
  a = rand(3,2)
  b = 3
  cache = kernel_cache(f,a,b)
  @time repeat(n,apply_kernel!,cache,f,a,b)
end

function bench4(n)
  f = elem(+)
  a = 5
  b = 3
  cache = kernel_cache(f,a,b)
  @time repeat(n,apply_kernel!,cache,f,a,b)
end

function bench5(n)
  f = elem(+)
  a = 3
  b = rand(3,2)
  cache = kernel_cache(f,a,b)
  @time repeat(n,apply_kernel!,cache,f,a,b)
end

function bench6(n)
  f = elem(+)
  a = rand(3,2)
  b = rand(3,2)
  cache = kernel_cache(f,a,b)
  @time repeat(n,apply_kernel!,cache,f,a,b)
end

function bench7(n)
  f = contract(-)
  a = rand(3,2)
  b = rand(3,2)
  cache = kernel_cache(f,a,b)
  @time repeat(n,apply_kernel!,cache,f,a,b)
end

for n in (1,1,10,1000,100000)
  @eval begin
    println("+++ runing suite for n = $($n) +++")
    bench1($n)
    bench2($n)
    bench3($n)
    bench4($n)
    bench5($n)
    bench6($n)
    bench7($n)
  end
end

end # module
