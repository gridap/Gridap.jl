# Gridap

[![Build Status](https://travis-ci.com/gridap/Gridap.jl.svg?branch=master)](https://travis-ci.com/gridap/Gridap.jl)
[![Codecov](https://codecov.io/gh/gridap/Gridap.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/Gridap.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/Gridap.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/Gridap.jl/dev) [![Join the chat at https://gitter.im/Gridap-jl/community](https://badges.gitter.im/Gridap-jl/community.svg)](https://gitter.im/Gridap-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## Installation 
```julia
# At least julia 1.1 required
using Pkg
Pkg.add("Gridap")
```
## Quick start

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
quad = CellQuadrature(trian,degree=2)

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
