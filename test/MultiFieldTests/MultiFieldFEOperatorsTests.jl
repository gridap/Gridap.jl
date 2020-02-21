module MultiFieldFEOperatorsTests

using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using SparseArrays
using Test
using LinearAlgebra

using Gridap.MultiField

order = 2

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)

strian = SkeletonTriangulation(model)
sdegree = order
squad = CellQuadrature(strian,sdegree)

V = TestFESpace(model=model,order=order,reffe=:Lagrangian,conformity=:H1,valuetype=Float64)
Q = TestFESpace(model=model,order=order-1,reffe=:Lagrangian,conformity=:L2,valuetype=Float64)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

function a(y,x)
  v,q = y
  u,p = x
  v*u + v*p - q*p
end

function l(y)
  v,q = y
  v*4 + q
end

function a_Γ(y,x)
  v,q = y
  u,p = x
  jump(v)*mean(u) + jump(∇(q))*jump(∇(p)) - mean(v)*mean(p)
end

t_Ω = AffineFETerm(a,l,trian,quad)
t_Γ = LinearFETerm(a_Γ,strian,squad)

op = AffineFEOperator(Y,X,t_Ω,t_Γ)

op = AffineFEOperator(SparseMatrixCSR,Y,X,t_Ω,t_Γ)

op = FEOperator(Y,X,t_Ω,t_Γ)
xh = zero(X)
b = residual(op,xh)
A = jacobian(op,xh)
test_fe_operator(op,get_free_values(xh),b)

op = FEOperator(SparseMatrixCSR,Y,X,t_Ω,t_Γ)

q = get_coordinates(quad)
w_q = get_weights(quad)
ϕ = get_cell_map(trian)
jac = ∇(ϕ)
jac_q = evaluate(jac,q)
x_q = evaluate(ϕ,q)

function cell_kernel!(A,B,y,x,j,w)

  A_vu = A[1,1]
  A_vp = A[1,2]
  A_qp = A[2,2]
  B_v = B[1]
  B_q = B[2]
  u = x[1]
  p = x[2]
  v = y[1]
  q = y[2]

  N_v, N_u = size(A_vu)
  N_q, N_p = size(A_qp)
  S = length(w)

  for s in 1:S
    dV = det(j[s])*w[s]
    for n_v in 1:N_v
      for n_u in 1:N_u
        A_vu[n_v,n_u] += v[s,n_v]*u[s,n_u]*dV
      end
      for n_p in 1:N_p
        A_vp[n_v,n_p] += v[s,n_v]*p[s,n_p]*dV
      end
    end
    for n_q in 1:N_q
      for n_p in 1:N_p
        A_qp[n_q,n_p] -=  q[s,n_q]*p[s,n_p]*dV
      end
    end
    for n_v in 1:N_v
      B_v[n_v] += v[s,n_v]*4*dV
    end
    for n_q in 1:N_q
      B_q[n_q] += q[s,n_q]*dV
    end
  end

end

function cellmat_Ω(y,x)
  y_q = evaluate(y,q)
  x_q = evaluate(x,q)
  apply_cellmatvec(cell_kernel!,y_q,x_q,jac_q,w_q)
end

t2_Ω = AffineFETermFromCellMatVec(cellmat_Ω,trian)

op = AffineFEOperator(Y,X,t_Ω)
op2 = AffineFEOperator(Y,X,t2_Ω)

@test get_matrix(op) == get_matrix(op2)
@test get_vector(op) == get_vector(op2)

end # module
