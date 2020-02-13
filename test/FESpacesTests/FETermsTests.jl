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

assem = SparseMatrixAssembler(V,U)

uhd = zero(U)
v = get_cell_basis(V)
u = get_cell_basis(U)
uh = interpolate(U,u_sol)

a(v,u) = ∇(v)*∇(u)
l(v) = f*v

j(u,v,du) = a(v,du)
r(u,v) = a(v,u) - l(v)

z(v,u) = jump(v)*jump(u)

t_linear = LinearFETerm(z,strian,squad)

w(v) = jump(v)*f
t_source = FESource(w,strian,squad)


# AffineFETerm

t_affine = AffineFETerm(a,l,trian,quad)

matdata = collect_cell_matrix(v,u,[t_affine, t_linear])
vecdata = collect_cell_vector(v,uhd,[t_affine, t_linear, t_source])
A = assemble_matrix(assem,matdata...)
b = assemble_vector(assem,vecdata...)
x = A \ b
@test x ≈ get_free_values(uh)

data = collect_cell_matrix_and_vector(v,u,uhd,[t_affine, t_linear, t_source])
A, b = allocate_matrix_and_vector(assem,data...)
assemble_matrix_and_vector!(A,b,assem,data...)
x = A \ b
@test x ≈ get_free_values(uh)
A, b = assemble_matrix_and_vector(assem,data...)
x = A \ b
@test x ≈ get_free_values(uh)

matdata = collect_cell_jacobian(uh,v,u,[t_affine])
vecdata = collect_cell_residual(uh,v,[t_affine])
A = assemble_matrix(assem,matdata...)
b = assemble_vector(assem,vecdata...)
x = A \ -b
@test (x.+1) ≈ ones(length(x))

# LinearFETerm and FESource

t_linear = LinearFETerm(a,trian,quad)
t_source = FESource(l,trian,quad)

matdata = collect_cell_matrix(v,u,[t_linear, t_source])
vecdata = collect_cell_vector(v,uhd,[t_linear, t_source])
A = assemble_matrix(assem,matdata...)
b = assemble_vector(assem,vecdata...)
x = A \ b
@test x ≈ get_free_values(uh)

matdata = collect_cell_jacobian(uh,v,u,[t_linear, t_source])
vecdata = collect_cell_residual(uh,v,[t_linear, t_source])
A = assemble_matrix(assem,matdata...)
b = assemble_vector(assem,vecdata...)
x = A \ -b
@test (x.+1) ≈ ones(length(x))

# NonLinearFETerm

t_nonlinear = FETerm(r,j,trian,quad)

matdata = collect_cell_jacobian(uh,v,u,[t_nonlinear])
vecdata = collect_cell_residual(uh,v,[t_nonlinear])
A = assemble_matrix(assem,matdata...)
b = assemble_vector(assem,vecdata...)
x = A \ -b
@test (x.+1) ≈ ones(length(x))

end # module
