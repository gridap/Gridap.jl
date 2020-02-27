module AffineFEOperatorsTests

using Test
using Gridap.Arrays
using Gridap.Algebra
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
using Gridap.Fields
using Gridap.FESpaces
using LinearAlgebra

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)
trian = get_triangulation(model)
degree = 4
quad = CellQuadrature(trian,degree)

order = 2
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

dirichlet_tags = [1,10]
V = GradConformingFESpace(reffes,model,dirichlet_tags)

U = TrialFESpace(V)

f(x) = x[2]

v = get_cell_basis(V)
u = get_cell_basis(U)

cellmat = integrate(∇(v)*∇(u),trian,quad)
cellvec = integrate(v*f,trian,quad)
cellids = collect(1:num_cells(trian))

assem = SparseMatrixAssembler(V,U)
A =  assemble_matrix(assem,[cellmat],[cellids],[cellids])
b =  assemble_vector(assem,[cellvec],[cellids])

op = AffineFEOperator(U,V,A,b)
@test A === get_matrix(op)
@test b === get_vector(op)

x = ones(length(b))
r = A*x - b

test_fe_operator(op,x,r,≈,jac=A)
@test isa(get_algebraic_operator(op), AffineOperator)

#

tol = 1.0e-9

u_sol(x) = x[1]+x[2]
f_fun(x) = 0

dirichlet_tags = "boundary"
V = GradConformingFESpace(reffes,model,dirichlet_tags)
U = TrialFESpace(V,u_sol)

a(v,u) = ∇(v)*∇(u)
l(v) = v*f_fun

t_Ω = AffineFETerm(a,l,trian,quad)

assem = SparseMatrixAssembler(V,U)

op = AffineFEOperator(U,V,assem,t_Ω)
uh = solve(op)
e = u_sol - uh
@test sum(integrate(e*e,trian,quad)) < tol

op = AffineFEOperator(U,V,t_Ω)
uh = solve(op)
e = u_sol - uh
@test sum(integrate(e*e,trian,quad)) < tol

q = get_coordinates(quad)
w_q = get_weights(quad)
ϕ = get_cell_map(trian)
jac = ∇(ϕ)
jac_q = evaluate(jac,q)
x_q = evaluate(ϕ,q)

function poisson_matvec_kernel!(mat,vec,∇v,∇u,v,j,w,x)
  Q = length(w)
  M,N = size(mat)
  for q in 1:Q
    dV = det(j[q])*w[q]
    f_q = f_fun(x[q])
    for n in 1:N
      for m in 1:M
        mat[m,n] += ∇v[q,m]*∇u[q,n]*dV
      end
    end
    for m in 1:M
      vec[m] += v[q,m]*f_q*dV
    end
  end
end

function cellmatvec_Ω(v,u)
  v_q = evaluate(v,q)
  ∇v_q = evaluate(∇(v),q)
  apply_cellmatvec(poisson_matvec_kernel!, ∇v_q, ∇v_q, v_q, jac_q, w_q, x_q)
end

t_Ω = AffineFETermFromCellMatVec(cellmatvec_Ω,trian)

op = AffineFEOperator(U,V,assem,t_Ω)
uh = solve(op)
e = u_sol - uh
@test sum(integrate(e*e,trian,quad)) < tol


end # module
