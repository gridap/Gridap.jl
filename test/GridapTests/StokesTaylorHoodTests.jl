module StokesTaylorHoodTests

using Test
using Gridap
import Gridap: ∇

#using LinearAlgebra: tr, ⋅

# Using automatic differentiation
u(x) = VectorValue( x[1]^2 + 2*x[2]^2, -x[1]^2 )
p(x) = x[1] + 3*x[2]
f(x) = -Δ(u)(x) + ∇(p)(x)
g(x) = (∇⋅u)(x)
∇u(x) = ∇(u)(x)

#u(x) = VectorValue( x[1]^2 + 2*x[2]^2, -x[1]^2 )
#∇u(x) = TensorValue( 2*x[1], 4*x[2], -2*x[1], zero(x[1]) )
#Δu(x) = VectorValue( 6, -2 )
#
#p(x) = x[1] + 3*x[2]
#∇p(x) = VectorValue(1,3)
#
#f(x) = -Δu(x) + ∇p(x)
#g(x) = tr(∇u(x))
#
#∇(::typeof(u)) = ∇u
#∇(::typeof(p)) = ∇p

domain = (0,2,0,2)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet",[1,2,5])
add_tag_from_tags!(labels,"neumann",[6,7,8])

order = 2

reffe_u = ReferenceFE(:Lagrangian,VectorValue{2,Float64},order)
reffe_p = ReferenceFE(:Lagrangian,Float64,order-1)

V = TestFESpace(model,reffe_u,labels=labels,dirichlet_tags="dirichlet",conformity=:H1)
Q = TestFESpace(model,reffe_p,labels=labels,conformity=:H1)

U = TrialFESpace(V,u)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model,labels,tags="neumann")
n_Γ = get_normal_vector(Γ)

degree = order
dΩ = LebesgueMeasure(Ω,degree)
dΓ = LebesgueMeasure(Γ,degree)

a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )*dΩ

l((v,q)) = ∫( v⋅f + q*g )*dΩ + ∫( v⋅(n_Γ⋅∇u) - (n_Γ⋅v)*p )*dΓ

op = AffineFEOperator(a,l,X,Y)

uh, ph = solve(op)

eu = u - uh
ep = p - ph

l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))

eu_l2 = l2(eu)
eu_h1 = h1(eu)
ep_l2 = l2(ep)

tol = 1.0e-9
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
