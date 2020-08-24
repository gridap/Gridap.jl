module FESolversTests

using Test
using Gridap.Arrays
using Gridap.Algebra
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
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
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

dirichlet_tags = [1,10]
V = GradConformingFESpace(reffes,model,dirichlet_tags)

U = TrialFESpace(V)

f(x) = x[2]

v = get_cell_basis(V)
u = get_cell_basis(U)

cellmat = integrate(∇(v)⊙∇(u),trian,quad)
cellvec = integrate(v⊙f,trian,quad)
cellids = collect(1:num_cells(trian))

assem = SparseMatrixAssembler(U,V)
matdata = ([cellmat],[cellids],[cellids])
vecdata = ([cellvec],[cellids])
A =  assemble_matrix(assem,matdata)
b =  assemble_vector(assem,vecdata)
x = A \ b
x0 = zeros(length(x))

op = AffineFEOperator(U,V,A,b)
solver = LinearFESolver()
test_fe_solver(solver,op,x0,x)
uh = solve(solver,op)
@test get_free_values(uh) ≈ x
uh = solve(op)
@test get_free_values(uh) ≈ x

solver = NonlinearFESolver()
test_fe_solver(solver,op,x0,x)
uh = solve(solver,op)
@test get_free_values(uh) ≈ x

end # module
