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
using Gridap.CellData

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

cellmat = integrate(∇(v)⊙∇(u),quad)
cellvec = integrate(v*f,quad)
cellids = collect(1:num_cells(trian))

assem = SparseMatrixAssembler(U,V)
matdata = ([cellmat],[cellids],[cellids])
vecdata = ([cellvec],[cellids])
A =  assemble_matrix(assem,matdata)
b =  assemble_vector(assem,vecdata)

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

a(v,u) = ∇(v)⊙∇(u)
l(v) = v*f_fun

t_Ω = AffineFETerm(a,l,trian,quad)

assem = SparseMatrixAssembler(V,U)

op = AffineFEOperator(U,V,assem,t_Ω)
uh = solve(op)
e = u_sol - uh
@test sum(integrate(e*e,quad)) < tol

op = AffineFEOperator(U,V,t_Ω)
uh = solve(op)
e = u_sol - uh
@test sum(integrate(e*e,quad)) < tol

end # module
