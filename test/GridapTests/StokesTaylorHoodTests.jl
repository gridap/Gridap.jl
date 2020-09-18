module StokesTaylorHoodTests

using Test
using Gridap
import Gridap: ∇

using LinearAlgebra: tr, ⋅

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

ref_style = [:reference,:physical]

for ref_st in ref_style
  V = TestFESpace(
  model=model,
  order=order,
  reffe=:Lagrangian,
  labels=labels,
  valuetype=VectorValue{2,Float64},
  dirichlet_tags="dirichlet",
  dof_space=ref_st,
  conformity=:H1)

  Q = TestFESpace(
  model=model,
  order=order-1,
  reffe=:Lagrangian,
  valuetype=Float64,
  dof_space=ref_st,
  conformity=:H1)

  U = TrialFESpace(V,u)
  P = TrialFESpace(Q)

  trian = get_triangulation(model)
  degree = order
  dΩ = LebesgueMeasure(trian,degree)

  btrian = BoundaryTriangulation(model,labels,"neumann")
  bdegree = order
  dΓ = LebesgueMeasure(btrian,bdegree)
  n = get_normal_vector(btrian)

  @form a((u,p),(v,q)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )*dΩ
  @form l((v,q)) = ∫( v⋅f + q*g )*dΩ + ∫( v⋅(n⋅∇u) - (n⋅v)*p )*dΓ

  uh, ph = solve( a((U,P),(V,Q))==l((V,Q)) )

  eu = u - uh
  ep = p - ph

  l2(v) = v⋅v
  h1(v) = v⋅v + ∇(v)⊙∇(v)

  eu_l2 = sqrt(sum(∫(l2(eu))*dΩ))
  eu_h1 = sqrt(sum(∫(h1(eu))*dΩ))
  ep_l2 = sqrt(sum(∫(l2(ep))*dΩ))

  tol = 1.0e-9
  @test eu_l2 < tol
  @test eu_h1 < tol
  @test ep_l2 < tol
end


end # module
