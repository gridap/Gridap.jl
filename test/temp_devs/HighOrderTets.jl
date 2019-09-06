module HighOrderTets

##
using Gridap, Test
using Gridap.CellValuesGallery

import Gridap: ∇

using Gridap.Helpers
using UnstructuredGrids.Kernels: refine_grid_connectivity
using UnstructuredGrids.Kernels: generate_tface_to_face
using Gridap.DiscreteModels: DiscreteModelFromData
##

##
ufun(x) = x[1] + x[2]
ufun_grad(x) = VectorValue(1.0,1.0,0.0)
∇(::typeof(ufun)) = ufun_grad
bfun(x) = 0.0

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0,0.0,1.0), partition=(4,4,4))
model = simplexify(model)


order = 4
# diritag = [1,2,3,4]
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

# Define test and trial
V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

# Define integration mesh and quadrature
trian = Triangulation(model)
quad = CellQuadrature(trian,order=6)

# Define forms
a(v,u) = inner(∇(v), ∇(u))
b(v) = inner(v,bfun)

uh = interpolate(U,ufun)
sum(integrate(a(uh,uh),trian,quad))

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

@test el2 < 1.e-8
@test eh1 < 1.e-8
##
uhno = uh
uho = uh
writevtk(trian,"trian",cellfields=["uh"=>uh])

# 1) Clean constructors without D or T

# Create dofbasis using node array for Lagrangian FEs

# Create BasisWithChangeOfBasis
# i.e., CanonicalBasis given DOFs

# nfacetoowndofs

# D = 1
#


# Closure n-face nodes
# Method that given the set of nodes and the nfacedofs, returns the
# nface dofs on the closure of the n-face (only sense for Lagrangian FEs)
##

end # module
