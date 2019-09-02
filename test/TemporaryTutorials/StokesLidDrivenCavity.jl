module StokesLidDrivenCavity

##
using Test
using Gridap
import Gridap: ∇

##
# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(100,100))

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
# fespace1 = ConformingFESpace(VectorValue{2,Float64},model,order,diritags1)
# fespace2 = ConformingFESpace(Float64,model,order-1,diritags2)
fespace1 = CLagrangianFESpace(VectorValue{2,Float64},model,order,diritags1,dirimasks1)
fespace2 = CLagrangianFESpace(Float64,model,order-1,diritags2,dirimasks2)
# Define test and trial
V1 = TestFESpace(fespace1)
V2 = TestFESpace(fespace2)
V = [V1, V2]

uD_1(x) = VectorValue(0.0,0.0)
uD_2(x) = VectorValue(1.0,0.0)
uD = [ (i == 6) ? uD_2 : uD_1 for i = 1:8 ]
U1 = TrialFESpace(fespace1,uD)
uD_p(x) = 0.0
U2 = TrialFESpace(fespace2,[uD_p])
U = [U1, U2]

# Define integration mesh and quadrature for volume
trian = Triangulation(model)
quad = CellQuadrature(trian,order=4)

divfun(x,∇u) = tr((∇u))
div(u) = CellBasis(trian,divfun,∇(u))


# Terms in the volume
a(v,u) = inner(∇(v[1]),∇(u[1])) - inner(div(v[1]),u[2]) + inner(v[2],div(u[1]))
t_Ω = LinearFETerm(a,trian,quad)
# assem = SparseMatrixAssembler(V,U)
# op = LinearFEOperator(V,U,assem,t_Ω)
# @fverdugo : Do you understand with there is an error when using LinearFEOperator?

# To circumvent the previous error I create a zero forcing term
u2fun(x) = 0.0
b2fun(x) = u2fun(x)
b2field = CellField(trian,b2fun)
b(v) =  inner(v[2],b2field)
t_Ω = AffineFETerm(a,b,trian,quad)
assem = SparseMatrixAssembler(V,U)
op = LinearFEOperator(V,U,assem,t_Ω)

# Define the FESolver
ls = LUSolver()
solver = LinearFESolver(ls)

# Solve!
uh = solve(solver,op)
writevtk(trian,"results",cellfields=["uh"=>uh[1],"ph"=>uh[2]])
##
end # module
