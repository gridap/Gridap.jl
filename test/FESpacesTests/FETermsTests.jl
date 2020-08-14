module FETermsTests

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

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 1
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

u_sol(x) = x[1]+x[2]
f(x) = 0

dirichlet_tags = "boundary"
V = GradConformingFESpace(reffes,model,dirichlet_tags)
U = TrialFESpace(V,u_sol)

trian = get_triangulation(model)
degree = 2
quad = CellQuadrature(trian,degree)

strian = SkeletonTriangulation(model)
degree = 2
squad = CellQuadrature(strian,degree)

assem = SparseMatrixAssembler(U,V)

uhd = zero(U)
v = get_cell_basis(V)
u = get_cell_basis(U)
uh = interpolate(u_sol,U)

a(u,v) = ∇(v)⊙∇(u)
l(v) = f*v

j(u,du,v) = a(du,v)
r(u,v) = a(u,v) - l(v)

z(u,v) = jump(v)*jump(u)
w(v) = jump(v)*f

t_linear = LinearFETerm(z,strian,squad)

t_source = FESource(w,strian,squad)

# AffineFETerm

t_affine = AffineFETerm(a,l,trian,quad)

matdata = collect_cell_matrix(u,v,[t_affine, t_linear])
vecdata = collect_cell_vector(uhd,v,[t_affine, t_linear, t_source])
A = assemble_matrix(assem,matdata)
b = assemble_vector(assem,vecdata)
x = A \ b
@test x ≈ get_free_values(uh)

data = collect_cell_matrix_and_vector(uhd,u,v,[t_affine, t_linear, t_source])
A, b = allocate_matrix_and_vector(assem,data)
assemble_matrix_and_vector!(A,b,assem,data)
x = A \ b
@test x ≈ get_free_values(uh)
A, b = assemble_matrix_and_vector(assem,data)
x = A \ b
@test x ≈ get_free_values(uh)

matdata = collect_cell_jacobian(uh,u,v,[t_affine])
vecdata = collect_cell_residual(uh,v,[t_affine])
A = assemble_matrix(assem,matdata)
b = assemble_vector(assem,vecdata)
x = A \ -b
@test (x.+1) ≈ ones(length(x))

# LinearFETerm and FESource

t_linear = LinearFETerm(a,trian,quad)
t_source = FESource(l,trian,quad)

matdata = collect_cell_matrix(u,v,[t_linear, t_source])
vecdata = collect_cell_vector(uhd,v,[t_linear, t_source])
A = assemble_matrix(assem,matdata)
b = assemble_vector(assem,vecdata)
x = A \ b
@test x ≈ get_free_values(uh)

matdata = collect_cell_jacobian(uh,u,v,[t_linear, t_source])
vecdata = collect_cell_residual(uh,v,[t_linear, t_source])
A = assemble_matrix(assem,matdata)
b = assemble_vector(assem,vecdata)
x = A \ -b
@test (x.+1) ≈ ones(length(x))

# NonlinearFETerm

t_nonlinear = FETerm(r,j,trian,quad)

matdata = collect_cell_jacobian(uh,u,v,[t_nonlinear])
vecdata = collect_cell_residual(uh,v,[t_nonlinear])
A = assemble_matrix(assem,matdata)
b = assemble_vector(assem,vecdata)
x = A \ -b
@test (x.+1) ≈ ones(length(x))

data = collect_cell_jacobian_and_residual(uh,u,v,[t_nonlinear])
A, b = allocate_matrix_and_vector(assem,data)
assemble_matrix_and_vector!(A,b,assem,data)
x = A \ -b
@test (x.+1) ≈ ones(length(x))

# AffineFETermFromCellMatVec

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
    f_q = f(x[q])
    for n in 1:N
      for m in 1:M
        mat[m,n] += ∇v[q,m]⊙∇u[q,n]*dV
      end
    end
    for m in 1:M
      vec[m] += v[q,m]*f_q*dV
    end
  end
end

function matvecfun(u,v)
  v_q = evaluate(v,q)
  ∇v_q = evaluate(∇(v),q)
  apply_cellmatvec(poisson_matvec_kernel!, ∇v_q, ∇v_q, v_q, jac_q, w_q, x_q)
end

t_matvec_Ω = AffineFETermFromCellMatVec(matvecfun,trian)

matdata = collect_cell_matrix(u,v,[t_matvec_Ω,])
vecdata = collect_cell_vector(uhd,v,[t_matvec_Ω,])
A = assemble_matrix(assem,matdata)
b = assemble_vector(assem,vecdata)
x = A \ b
@test x ≈ get_free_values(uh)

data = collect_cell_matrix_and_vector(uhd,u,v,[t_matvec_Ω,])
A, b = allocate_matrix_and_vector(assem,data)
assemble_matrix_and_vector!(A,b,assem,data)
x = A \ b
@test x ≈ get_free_values(uh)

function poisson_jacres_kernel!(jac,res,∇v,∇du,v,∇uh,j,w,x)
  Q = length(w)
  M,N = size(jac)
  for q in 1:Q
    dV = det(j[q])*w[q]
    f_q = f(x[q])
    for m in 1:M
      for n in 1:N
        jac[m,n] += ∇v[q,m]⊙∇du[q,n]*dV
      end
      res[m] += ∇v[q,m]⊙∇uh[q]*dV
    end
    for m in 1:M
      res[m] -= v[q,m]*f_q*dV
    end
  end
end

function jacresfun(uh,du,v)
  v_q = evaluate(v,q)
  ∇v_q = evaluate(∇(v),q)
  ∇du_q = ∇v_q
  ∇uh_q = evaluate(∇(uh),q)
  apply_cellmatvec(poisson_jacres_kernel!, ∇v_q, ∇du_q, v_q, ∇uh_q, jac_q, w_q, x_q)
end

t_jacres_Ω = FETermFromCellJacRes(jacresfun,trian)

matdata = collect_cell_jacobian(uh,u,v,[t_jacres_Ω])
vecdata = collect_cell_residual(uh,v,[t_jacres_Ω])
A = assemble_matrix(assem,matdata)
b = assemble_vector(assem,vecdata)
x = A \ -b
@test (x.+1) ≈ ones(length(x))

data = collect_cell_jacobian_and_residual(uh,u,v,[t_jacres_Ω])
A, b = allocate_matrix_and_vector(assem,data)
assemble_matrix_and_vector!(A,b,assem,data)
x = A \ -b
@test (x.+1) ≈ ones(length(x))

end # module
