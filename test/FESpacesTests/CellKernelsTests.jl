module CellKernelsTests

using Test
using Gridap.Fields
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Integration
using Gridap.FESpaces
using LinearAlgebra
using Gridap.TensorValues

function poisson_matvec_kernel!(mat,vec,∇u,∇v,v,j,w)
  Q = length(w)
  M,N = size(mat)
  for q in 1:Q

    dV = det(j[q])*w[q]

    for n in 1:N
      for m in 1:M
        mat[m,n] += ∇v[q,m]⊙∇u[q,n]*dV
      end
    end

    for m in 1:M
      vec[m] += v[q,m]*dV
    end

  end
end

function poisson_mat_kernel!(mat,∇u,∇v,j,w)
  Q = length(w)
  M,N = size(mat)
  for q in 1:Q

    dV = det(j[q])*w[q]

    for n in 1:N
      for m in 1:M
        mat[m,n] += ∇v[q,m]⊙∇u[q,n]*dV
      end
    end

  end
end

function poisson_vec_kernel!(vec,v,j,w)
  Q = length(w)
  M = length(vec)
  for q in 1:Q

    dV = det(j[q])*w[q]

    for m in 1:M
      vec[m] += v[q,m]*dV
    end

  end
end

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 1
V = TestFESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=order,conformity=:H1)
U = TrialFESpace(V)

trian = Triangulation(model)
quad = CellQuadrature(trian,order)

v = get_cell_basis(V)
u = get_cell_basis(U)

q = get_coordinates(quad)
w = get_weights(quad)
j = ∇(get_cell_map(trian))

v_q = evaluate(v,q)
∇v_q = evaluate(∇(v),q)
j_q = evaluate(j,q)

cellmatvec = apply_cellmatvec(poisson_matvec_kernel!, ∇v_q, ∇v_q, v_q, j_q, w)
cellmat = apply_cellmatrix(poisson_mat_kernel!, ∇v_q, ∇v_q, j_q, w)
cellvec = apply_cellvector(poisson_vec_kernel!, v_q, j_q, w)

a(v,u) = ∇(v)⊙∇(u)
l(v) = v

cellmat2 = integrate(a(v,u),trian,quad)
cellvec2 = integrate(l(v),trian,quad)
cellmatvec2 = pair_arrays(cellmat2,cellvec2)

@test cellmatvec == cellmatvec2
@test cellmat == cellmat2
@test cellvec == cellvec2

end # module
