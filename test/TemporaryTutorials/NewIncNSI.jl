module PlayingWithFEMachinery

##
using Test
using Gridap
import Gridap: ∇

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# order = 2
# T = VectorValue{2,Float64}
# diritags = [1,2,5]
# dirimasks = [(true,false),(false,true),(true,true)]
# fespace = CLagrangianFESpace(T,model,order,diritags,dirimasks)
#
# T = Float64
# diritags = [1,2,5]
# dirimasks = [(true,),(false,),(true,)]
# fespace = CLagrangianFESpace(T,model,order,diritags,dirimasks)
#
#
# num_comps = 1
# dirimasks = Vector{NTuple{num_comps,Bool}}(undef,length(diritags))
# flag = ntuple(i -> true, num_comps)
# fill!(dirimasks,flag)
# fespace = CLagrangianFESpace(T,model,order,diritags,dirimasks)
#
# num_comps = 2
# dirimasks = Vector{NTuple{num_comps,Bool}}(undef,length(diritags))
# flag = ntuple(i -> true, num_comps)
# fill!(dirimasks,flag)
# fespace = CLagrangianFESpace(T,model,order,diritags,dirimasks)
#
# # with diritags it does not work
# # diritags = ["boundary"]
# diritags = [1]
# num_comps = 2
# dirimasks = Vector{NTuple{num_comps,Bool}}(undef,length(diritags))
# flag = ntuple(i -> true, num_comps)
# fill!(dirimasks,flag)
# fespace = CLagrangianFESpace(T,model,order,diritags,dirimasks)
#
# # if no boundary conditions are to be defined
# diritags = Int[]
# dirimasks = []
# fespace = DLagrangianFESpace(T,model,order,diritags,dirimasks)

# Now I define two spaces, one vector space of order 2
# and one scalar space of order 1
order = 2
T = VectorValue{2,Float64}
diritags = [1,2,3,4,5,6,7,8]
dirimasks = Vector{NTuple{2,Bool}}(undef,length(diritags))
flag = (true, true)
fill!(dirimasks,flag); dirimasks[8] = (true, false)
_Uh = ConformingFESpace(Float64,model,order,diritags)
# _Uh = CLagrangianFESpace(T,model,order,diritags)
uD_1(x) = VectorValue(0.0,0.0)
uD_2(x) = VectorValue(1.0,0.0)
uD = [ (i == 8) ? uD_2 : uD_1 for i = 1:8 ]

Uh = TrialFESpace(_Uh,uD)
Vh = TestFESpace(_Uh)

order = 1
T = Float64
diritags = [1,]
dirimasks = [(true,)]
_Ph = DLagrangianFESpace(T,model,order,diritags,dirimasks)
pD(x) = 0.0; pD_grad(x) = VectorValue(0.0,0.0)
Ph = TrialFESpace(_Ph,pD)
Qh = TestFESpace(_Ph)

Xh = [Uh, Ph]
Yh = [Vh, Qh]

a(u,v) = inner(∇(v[1]),∇(u[1])) + inner(v[2],u[2])
da(u,v,du) = a(du,v)
res(u,v) = a(u,v)
jac(u,v,du) = da(u,v,du)

# Define Assembler
assem = SparseMatrixAssembler(Yh,Xh)

# Define the FEOperator
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)
op = NonLinearFEOperator(res,jac,Yh,Xh,assem,trian,quad)

# Define the FESolver
ls = LUSolver()
tol = 1.e-10
maxiters = 20
nls = NewtonRaphsonSolver(ls,tol,maxiters)
solver = NonLinearFESolver(nls)

# Solve!
uh = solve(solver,op)
end
