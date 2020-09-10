module FEOperatorsFromTermsTests

using Test
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Algebra
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
using Gridap.Fields
using Gridap.FESpaces
using Gridap.CellData

import Gridap.Fields: ∇

u(x) = x[1] + x[2]
ufun_grad(x) = VectorValue(1.0,1.0)
∇(::typeof(u)) = ufun_grad
f(x) = -(3.0*x[1]+x[2]+1.0)

@law ν(u,x) = (u+1.0)*x[1]
@law dν(du,x) = x[1]*du

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

order = 2
diritag = "boundary"
V = GradConformingFESpace(Float64,model,order,diritag)

U = TrialFESpace(V,u)

trian = get_triangulation(model)
degree = 3*order-1
quad = CellQuadrature(trian,degree)

const x = get_physical_coordinate(trian) # its safer to mark as constant objects that are captured

a(u,v,du) = inner( ∇(v), ν(u,x)*∇(du))
res(u,v) = a(u,v,u) - inner(v,f)
jac(u,v,du) = a(u,v,du) + inner(∇(v),dν(du,x)*∇(u))

t_Ω = FETerm(res,jac,trian,quad)
op = FEOperator(U,V,t_Ω)

uh = solve(op)

e = u - uh

l2(u) = inner(u,u)
sh1(u) = inner(∇(u),∇(u))
h1(u) = sh1(u) + l2(u)

el2 = sqrt(sum( integrate(l2(e),quad) ))
eh1 = sqrt(sum( integrate(h1(e),quad) ))
ul2 = sqrt(sum( integrate(l2(uh),quad) ))
uh1 = sqrt(sum( integrate(h1(uh),quad) ))

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

y = ones(num_free_dofs(U))
uh = FEFunction(U,y)
r = residual(op,uh)
A = jacobian(op,uh)
test_fe_operator(op,y,r,≈,jac=A)

end # module
