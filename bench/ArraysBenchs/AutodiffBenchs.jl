module AutodiffBenchs

using Gridap.Arrays
import Gridap.Arrays: kernel_cache
import Gridap.Arrays: apply_kernel!

function user_cell_energy(cell_u)
  function f(u)
    e = zero(eltype(u))
    for ui in u
      e += ui^2
    end
    e
  end
  apply(f,cell_u)
end

function user_cell_residual(cell_u)
  apply(R(),cell_u)
end

struct R <: Kernel end

function kernel_cache(k::R,u)
  r = copy(u)
  r .= zero(eltype(u))
  r
end

@inline function apply_kernel!(r,k::R,u)
  for (i,ui) in enumerate(u)
    r[i] = 2*ui
  end
  r
end

function user_cell_jacobian(cell_u)
  apply(J(),cell_u)
end

struct J <: Kernel end

function kernel_cache(k::J,u)
  n = length(u)
  j = zeros(n,n)
  j
end

@inline function apply_kernel!(j,k::J,u)
  n = length(u)
  for i in 1:n
    j[i,i] = 2
  end
  j
end

L = 10^6
l = 14
cell_u = [ rand(l) for i in 1:L ]

cell_e = user_cell_energy(cell_u)
cell_r = user_cell_residual(cell_u)
cell_j = user_cell_jacobian(cell_u)
cell_h = cell_j

cell_r_auto = autodiff_array_gradient(user_cell_energy,cell_u)
cell_j_auto = autodiff_array_jacobian(user_cell_residual,cell_u)
cell_h_auto = autodiff_array_hessian(user_cell_energy,cell_u)

function bench_loop(a)
  cache = array_cache(a)
  @time loop_cached(a,cache)
end

function loop_cached(a,cache)
  for i in eachindex(a)
    ai = getindex!(cache,a,i)
  end
end

println("Warm up ...")
#bench_loop(cell_e)
bench_loop(cell_r)
bench_loop(cell_r_auto)
bench_loop(cell_j)
bench_loop(cell_j_auto)
bench_loop(cell_h)
bench_loop(cell_h_auto)


println("Run ...")
#bench_loop(cell_e)
bench_loop(cell_r)
bench_loop(cell_r_auto)
bench_loop(cell_j)
bench_loop(cell_j_auto)
bench_loop(cell_h)
bench_loop(cell_h_auto)

#using BenchmarkTools
#
#cache = array_cache(cell_h_auto)
#i = 3
#@btime a = getindex!($cache,$cell_h_auto,$i)

end # module
