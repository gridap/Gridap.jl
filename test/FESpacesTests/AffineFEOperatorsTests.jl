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

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(model,reffe;dirichlet_tags=[1,10])
U = V

trian = get_triangulation(model)
degree = 4
quad = CellQuadrature(trian,degree)

#f(x) = x[2]
#v = get_cell_shapefuns(V)
#u = get_cell_shapefuns_trial(U)
#
#cellmat = integrate(∇(v)⊙∇(u),quad)
#cellvec = integrate(v*f,quad)
#cellids = collect(1:num_cells(trian))
#
#assem = SparseMatrixAssembler(U,V)
#matdata = ([cellmat],[cellids],[cellids])
#vecdata = ([cellvec],[cellids])
#A =  assemble_matrix(assem,matdata)
#b =  assemble_vector(assem,vecdata)
#
#op = AffineFEOperator(U,V,A,b)
#@test A === get_matrix(op)
#@test b === get_vector(op)
#
#x = ones(length(b))
#r = A*x - b
#
#test_fe_operator(op,x,r,≈,jac=A)
#@test isa(get_algebraic_operator(op), AffineOperator)

#

tol = 1.0e-9

u_sol(x) = x[1]+x[2]
f_fun(x) = 0

V = FESpace(model,reffe;dirichlet_tags="boundary")
U = TrialFESpace(V,u_sol)

dΩ = Measure(quad)

a(u,v) = ∫(∇(v)⊙∇(u))*dΩ
l(v) = ∫(v*f_fun)*dΩ

assem = SparseMatrixAssembler(U,V)

op = AffineFEOperator(U,V,assem) do u,v
  ∫(∇(v)⊙∇(u))*dΩ, ∫(v*f_fun)*dΩ
end
uh = solve(op)
e = u_sol - uh

@test sum(∫(e*e)*dΩ) < tol
@test ∑(∫(e*e)*dΩ) < tol

op = AffineFEOperator(U,V) do u,v
  ∫(∇(v)⊙∇(u))*dΩ, ∫(v*f_fun)*dΩ
end
uh = solve(op)
e = u_sol - uh
@test sum(∫(e*e)*dΩ) < tol

op = AffineFEOperator(a,l,U,V)
uh = solve(op)
e = u_sol - uh
@test sum(∫(e*e)*dΩ) < tol

end # module
