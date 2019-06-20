# Getting Started

## Installation requirements

At the moment, Gridap requires at least Julia version 1.1.

Gridap has been tested on Linux and Mac Os operating systems.

## Installation

Gridap is a registered package. Thus, the installation should be straight forward
using the Julia's package manager [Pkg](https://julialang.github.io/Pkg.jl/v1/):

```julia
using Pkg
Pkg.add("Gridap")
```

For further information about how to install and manage Julia packages, see the
[Pkg documentation](https://julialang.github.io/Pkg.jl/v1/).

## First example

Solve a Poisson problem on the unit square with Dirichlet boundary conditions

```julia
using Gridap
import Gridap: ∇

# Define manufactured functions
ufun(x) = x[1] + x[2]
ufun_grad(x) = VectorValue(1.0,1.0)
∇(::typeof(ufun)) = ufun_grad
bfun(x) = 0.0

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

# Construct the FEspace
order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

# Define test and trial spaces
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Define integration mesh and quadrature
trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

# Define the source term
bfield = CellField(trian,bfun)

# Define forms of the problem
a(v,u) = inner(∇(v), ∇(u))
b(v) = inner(v,bfield)

# Define Assembler
assem = SparseMatrixAssembler(V,U)

# Define the FEOperator
op = LinearFEOperator(a,b,V,U,assem,trian,quad)

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

@assert el2 < 1.e-8
@assert eh1 < 1.e-8

# Write the numerical solution, the manufactured solution, and the error
# in a vtu file
writevtk(trian,"results",nref=4,cellfields=["uh"=>uh,"u"=>u,"e"=>e])
```

