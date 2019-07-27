module VectorValuedFEOperatorsTests

using Test
using Gridap
import Gridap: ∇

# Define manufactured functions
ufun(x) = VectorValue(x[1] + x[2],x[1])
ufun_grad(x) = TensorValue(1.0,1.0,1.0,0.0)
∇(::typeof(ufun)) = ufun_grad
bfun(x) = VectorValue(0.0,0.0)
gfun(x) = VectorValue(1.0,1.0)

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace
order = 2
diritags = [1,2,3,4,5,6,7]
T = VectorValue{2,Float64}
fespace = ConformingFESpace(T,model,order,diritags)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Setup integration on the volume
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

# Setup integration on Neumann boundary
neumanntags = [8,]
btrian = BoundaryTriangulation(model,neumanntags)
bquad = CellQuadrature(btrian,order=2)

# Terms in the volume
bfield = CellField(trian,bfun)
a(v,u) = inner(∇(v), ∇(u))
b(v) = inner(v,bfield)
t_Ω = AffineFETerm(a,b,trian,quad)

# Terms on the Neumann Boundary
gfield = CellField(btrian,gfun)
g(v) = inner(v,gfield)
t_ΓN = FESource(g,btrian,bquad)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = LinearFEOperator(V,U,assem,t_Ω,t_ΓN)

# Define the FESolver
ls = LUSolver()
solver = LinearFESolver(ls)

# Solve!
uh = solve(solver,op)

# Define exact solution and error
u = CellField(trian,ufun)
e = u - uh

# Define norms to measure the error
l2(u) = inner(u,u)
h1(u) = a(u,u) + l2(u)

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))

@test el2 < 1.e-8
@test eh1 < 1.e-8

# Idem for linear elasticity

# Material parameters
const E = 10
const ν = 0.25
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))

# Constitutive law
σfun(x,ε) = λ*tr(ε)*one(ε) + 2*μ*ε

# Define manufactured functions
ufun(x) = VectorValue(x[1] + x[2],x[1])
ufun_grad(x) = TensorValue(1.0,1.0,1.0,0.0)
∇(::typeof(ufun)) = ufun_grad
bfun(x) = VectorValue(0.0,0.0)

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace
order = 2
diritag = "boundary"
T = VectorValue{2,Float64}
fespace = ConformingFESpace(T,model,order,diritag)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Setup integration on the volume
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

# Terms in the volume
bfield = CellField(trian,bfun)
σ(ε) = CellBasis(trian,σfun,ε)
a_elast(v,u) = inner( ε(v), σ(ε(u)) )
b(v) = inner(v,bfield)
t_Ω = AffineFETerm(a_elast,b,trian,quad)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = LinearFEOperator(V,U,assem,t_Ω)

# Define the FESolver
ls = LUSolver()
solver = LinearFESolver(ls)

# Solve!
uh = solve(solver,op)

# Define exact solution and error
u = CellField(trian,ufun)
e = u - uh

# Compute errors
el2 = sqrt(sum( integrate(l2(e),trian,quad) ))
eh1 = sqrt(sum( integrate(h1(e),trian,quad) ))
eener = sqrt(sum( integrate(a_elast(e,e),trian,quad) ))

@test el2 < 1.e-8
@test eh1 < 1.e-8
@test eener < 1.e-8

end # module
