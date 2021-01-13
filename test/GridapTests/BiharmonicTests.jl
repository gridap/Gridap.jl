module BiharmonicTests

using Test
using Gridap

# Analytical manufactured solution
u(x) = x[1]*(x[1]-1)*x[2]*(x[2]-1)
f(x) = Δ(Δ(u))(x)
g(x) = Δ(u)(x)
@test f(VectorValue(0.5,0.5)) == 8.0
@test g(VectorValue(0.5,0.5)) == -1

# Domain
domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

# FE space
order = 2
V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order),dirichlet_tags="boundary")
U = TrialFESpace(V,u)

# Triangulation
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
dΛ = Measure(Λ,degree)
nΓ = get_normal_vector(Γ)
nΛ = get_normal_vector(Λ)

# Weak form
const h = (domain[2]-domain[1]) / partition[1]
const γ = 1
a(u,v) = ∫( Δ(u)*Δ(v) )dΩ +
         ∫( - mean(Δ(u))*jump(∇(v)⋅nΛ) - jump(∇(u)⋅nΛ)*mean(Δ(v)) + γ/h*jump(∇(u)⋅nΛ)*jump(∇(v)⋅nΛ) )dΛ
l(v) = ∫( v*f )dΩ + ∫( g*(∇(v)⋅nΓ) )dΓ
op = AffineFEOperator(a,l,U,V)

uₕ = solve(op)

# Error
e = u - uₕ
l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))
el2 = l2(e)
eh1 = h1(e)
tol = 1.0e-10
@test el2 < tol
@test eh1 < tol

end
