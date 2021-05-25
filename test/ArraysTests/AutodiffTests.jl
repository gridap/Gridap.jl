module AutodiffTests

using Test

using Gridap.Arrays

function user_cell_energy(cell_u)
  function f(u)
    e = zero(eltype(u))
    for ui in u
      e += ui^2
    end
    e
  end
  lazy_map(f,cell_u)
end

function user_cell_residual(cell_u)
  lazy_map(Broadcasting(i->2*i),cell_u)
end

function user_cell_jacobian(cell_u)
  function f(u)
    n = length(u)
    j = zeros(n,n)
    for i in 1:n
      j[i,i] = 2
    end
    j
  end
  lazy_map(f,cell_u)
end

L = 10
l = 4 # Do not use a number <13 (too easy for ForwardDiff)

cell_u = [rand(l) for i in 1:L]
cell_e = user_cell_energy(cell_u)
cell_r = user_cell_residual(cell_u)
cell_j = user_cell_jacobian(cell_u)
cell_h = cell_j

cell_j_auto = autodiff_array_jacobian(user_cell_residual,cell_u)
cell_r_auto = autodiff_array_gradient(user_cell_energy,cell_u)
cell_h_auto = autodiff_array_hessian(user_cell_energy,cell_u)
test_array(cell_r_auto,cell_r)
test_array(cell_j_auto,cell_j)
test_array(cell_h_auto,cell_h)

ids = [3,4,1,2]

function user_cell_energy_Γ(cell_u)
  cell_u_Γ = lazy_map(Reindex(cell_u),ids)
  user_cell_energy(cell_u_Γ)
end

function user_cell_residual_Γ(cell_u)
  cell_u_Γ = lazy_map(Reindex(cell_u),ids)
  user_cell_residual(cell_u_Γ)
end

function user_cell_jacobian_Γ(cell_u)
  cell_u_Γ = lazy_map(Reindex(cell_u),ids)
  user_cell_jacobian(cell_u_Γ)
end

cell_e_Γ = user_cell_energy_Γ(cell_u)
cell_r_Γ = user_cell_residual_Γ(cell_u)
cell_j_Γ = user_cell_jacobian_Γ(cell_u)
cell_h_Γ = cell_j_Γ

cell_r_Γ_auto = autodiff_array_gradient(user_cell_energy_Γ,cell_u,ids)
cell_j_Γ_auto = autodiff_array_jacobian(user_cell_residual_Γ,cell_u,ids)
cell_h_Γ_auto = autodiff_array_hessian(user_cell_energy_Γ,cell_u,ids)

test_array(cell_r_Γ_auto,cell_r_Γ)
test_array(cell_j_Γ_auto,cell_j_Γ)
test_array(cell_h_Γ_auto,cell_h_Γ)

end # module
