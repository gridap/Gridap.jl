module PLaplacianTests

using Test
using Gridap
import Gridap: ∇
using LinearAlgebra

domain = (0,1,0,1)
cells = (4,4)
model = CartesianDiscreteModel(domain,cells)

u(x) = x[1] + x[2]
∇u(x) = VectorValue(1,1)
f(x) = 0
∇(::typeof(u)) = ∇u

const p = 3
flux(∇u) = norm(∇u)^(p-2) * ∇u
dflux(∇du,∇u) = (p-2)*norm(∇u)^(p-4)*inner(∇u,∇du)*∇u + norm(∇u)^(p-2)*∇du

order = 3
T = Float64

V = TestFESpace(model,ReferenceFE(:Lagrangian,Float64,order),dirichlet_tags="boundary")
U = TrialFESpace(V,u)

Ω = Triangulation(model)

degree = 2*order
dΩ = LebesgueMeasure(Ω,degree)

res(u,v) = ∫( ∇(v)⋅(flux∘∇(u)) )*dΩ
jac(u,du,v) = ∫( ∇(v)⋅(dflux∘(∇(du),∇(u))) )*dΩ

op = FEOperator(res,jac,U,V)

nls = NLSolver(show_trace=false, method=:newton)
solver = FESolver(nls)

x = rand(T,num_free_dofs(U))
uh0 = FEFunction(U,x)
uh, = solve!(uh0,solver,op)

e = u - uh

l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

# now with autodiff

op = FEOperator(res,U,V)

x = rand(T,num_free_dofs(U))
uh0 = FEFunction(U,x)
uh, = solve!(uh0,solver,op)

e = u - uh

l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

end # module
