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

op = AffineFEOperator(V,U,A,b)
@test A === get_matrix(op)
@test b === get_vector(op)

x = ones(length(b))
r = A*x - b

test_fe_operator(op,x,r,≈,jac=A)

end # module
