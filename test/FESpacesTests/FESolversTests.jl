module FESolversTests

using Test
using Gridap.Arrays
using Gridap.Algebra
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Fields
using Gridap.FESpaces
using Gridap.CellData

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)
trian = get_triangulation(model)
degree = 4
quad = CellQuadrature(trian,degree)

order = 2

order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(model,reffe,dirichlet_tags=[1,10])
U = V
f(x) = x[2]

v = get_fe_basis(V)
u = get_trial_fe_basis(U)

cellmat = integrate(∇(v)⊙∇(u),quad)
cellvec = integrate(v⊙f,quad)
rows = get_cell_dof_ids(V,trian)
cols = get_cell_dof_ids(U,trian)
cellmat_c = attach_constraints_cols(U,cellmat,trian)
cellmat_rc = attach_constraints_rows(V,cellmat_c,trian)
cellvec_r = attach_constraints_rows(V,cellvec,trian)

assem = SparseMatrixAssembler(U,V)
matdata = ([cellmat_rc],[rows],[cols])
vecdata = ([cellvec_r],[rows])
A =  assemble_matrix(assem,matdata)
b =  assemble_vector(assem,vecdata)
x = A \ b
x0 = zeros(length(x))

op = AffineFEOperator(U,V,A,b)
solver = LinearFESolver()
test_fe_solver(solver,op,x0,x)
uh = solve(solver,op)
@test get_free_dof_values(uh) ≈ x
uh = solve(op)
@test get_free_dof_values(uh) ≈ x

solver = NonlinearFESolver()
test_fe_solver(solver,op,x0,x)
uh = solve(solver,op)
@test get_free_dof_values(uh) ≈ x

# Now using algebraic solvers directly
solver = LUSolver()
uh = solve(solver,op)
@test get_free_dof_values(uh) ≈ x
uh = solve(op)
@test get_free_dof_values(uh) ≈ x
uh,cache = solve!(uh,solver,op)
@test get_free_dof_values(uh) ≈ x
uh, = solve!(uh,solver,op,cache)

solver = NLSolver(LUSolver(),show_trace=false,method=:newton)
uh = solve(solver,op)
@test get_free_dof_values(uh) ≈ x
uh = solve(op)
@test get_free_dof_values(uh) ≈ x
uh,cache = solve!(uh,solver,op)
@test get_free_dof_values(uh) ≈ x
uh, = solve!(uh,solver,op,cache)


end # module
