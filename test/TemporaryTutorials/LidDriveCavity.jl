
module StokesTutorial

using Test
using Gridap
import Gridap: ∇

##
# Define manufactured functions
u1fun(x) = VectorValue(1.0,0.0)
u2fun(x) = VectorValue(0.0,0.0)
u2fun(x) = 0.0

u1fun_grad(x) = TensorValue(0.0,0.0,0.0,0.0)
u2fun_grad(x) = VectorValue(0.0,0.0)

∇(::typeof(u1fun)) = u1fun_grad
∇(::typeof(u2fun)) = u2fun_grad


b1fun(x) = u2fun_grad(x)
b2fun(x) = u2fun(x)

g1fun(x) = 1.0

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace.
# We allow each one of the fields to have different boundary conditions
order = 2
diritags1 = [1,2,3,4,5,6,7,8]
num_comps = 2
dirimasks1 = Vector{NTuple{num_comps,Bool}}(undef,length(diritags1))
flag = ntuple(i -> true, num_comps)
fill!(dirimasks1,flag)
diritags2 = [1]
dirimasks2 = [(true,)]
# This method is not stable but I can compute the matrix etc
fespace1 = ConformingFESpace(VectorValue{2,Float64},model,order,diritags1)
fespace2 = ConformingFESpace(Float64,model,order-1,diritags2)
# With these FE spaces there are bugs
fespace1 = CLagrangianFESpace(VectorValue{2,Float64},model,order,diritags1,dirimasks1)
fespace2 = DLagrangianFESpace(Float64,model,order-1,diritags2,dirimasks2)
##
# Define test and trial
V1 = TestFESpace(fespace1)
V2 = TestFESpace(fespace2)
V = [V1, V2]

uD_1(x) = VectorValue(0.0,0.0)
uD_2(x) = VectorValue(1.0,0.0)
uD = [ (i == 6) ? uD_2 : uD_1 for i = 1:8 ]
U1 = TrialFESpace(fespace1,uD)
U2 = TrialFESpace(fespace2,u2fun)
U = [U1, U2]
##

# Define integration mesh and quadrature for volume
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

b1field = CellField(trian,b1fun)
b2field = CellField(trian,b2fun)

divfun(x,∇u) = tr((∇u))
div(u) = CellBasis(trian,divfun,∇(u))


# Terms in the volume
a(v,u) = inner(∇(v[1]),∇(u[1])) + inner(v[2],u[2]) - inner(div(v[1]),u[2]) + inner(v[2],div(u[1]))
# a(u,v) = inner(∇(v[1]),ν*∇(u[1])) - inner(tr(∇(v[1])),u[2]) + inner(v[2],tr(∇(u[1])))
b(v) = inner(v[1],b1field) + inner(v[2],b2field)
t_Ω = AffineFETerm(a,b,trian,quad)

# Terms on Neumann boundary
# Note that the Neumann BC only applies on the first field
# g1field = CellField(btrian,g1fun)
# g(v) = inner(v[1],g1field)
# t_Γ = FESource(g,btrian,bquad)
##

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
# op = LinearFEOperator(V,U,assem,t_Ω,t_Γ)
op = LinearFEOperator(V,U,assem,t_Ω)

# Define the FESolver
ls = LUSolver()
solver = LinearFESolver(ls)
##

# Solve!
uh = solve(solver,op)
uh
# Define exact solution and error
u1 = CellField(trian,u1fun)
e1 = u1 - uh[1]
#
u2 = CellField(trian,u2fun)
e2 = u2 - uh[2]
#
# # Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = inner(∇(u),∇(u)) + l2(u)
#
# # Compute errorsConformingFESpace
e1l2 = sqrt(sum( integrate(l2(e1),trian,quad) ))
e1h1 = sqrt(sum( integrate(h1(e1),trian,quad) ))
#
e2l2 = sqrt(sum( integrate(l2(e2),trian,quad) ))
e2h1 = sqrt(sum( integrate(h1(e2),trian,quad) ))
##
# @test e1l2 < 1.e-8
# @test e1h1 < 1.e-8
# #
# @test e2l2 < 1.e-8
# @test e2h1 < 1.e-8

writevtk(trian,"results",nref=4,cellfields=["uh"=>uh[1],"ph"=>uh[2]])

end
##
#
#
# # Further tests
#
# # This is only to stress the single term API
# op = LinearFEOperator(a,b,V,U,assem,trian,quad)

end
