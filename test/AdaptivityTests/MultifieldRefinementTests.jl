module MultifieldGridTransferTests

using Test
using Gridap
using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays

f₁(x) = x[1]+x[2]
f₂(x) = x[1]

# Get refined model and triangulation
cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model = refine(cart_model,2)
trian = Triangulation(model)
ctrian = Triangulation(get_parent(model))

# Measures
dΩh   = Measure(trian,2)
dΩH   = Measure(ctrian,2)
dΩHh  = Measure(ctrian,trian,1)

# Coarse FESpace
et = Float64
reffe = LagrangianRefFE(et, QUAD, 1)
V₁    = FESpace(cart_model, reffe, conformity=:H1)
V₁²   = MultiFieldFESpace([V₁,V₁])
uH    = interpolate([f₁, f₂], V₁²)

# Fine FESpace
reffe = LagrangianRefFE(et, QUAD, 1)
V₂    = FESpace(model, reffe, conformity=:H1)
V₂²   = MultiFieldFESpace([V₂,V₂])
uh    = interpolate([f₁, f₂], V₂²)

# Interpolation
uHi = interpolate(uh, V₁²)
uhi = interpolate(uH, V₂²)

pts = [VectorValue(rand(2)) for i=1:10]
uHi₁,uHi₂ = uHi
uhi₁,uhi₂ = uhi
for pt in pts
  @test uhi₁(pt) ≈ f₁(pt)
  @test uhi₂(pt) ≈ f₂(pt)
  @test uHi₁(pt) ≈ f₁(pt)
  @test uHi₂(pt) ≈ f₂(pt)
end

# Integration
a((u,p),(v,q),dΩ) = ∫( v*u + p*u - q*p - v*4 + q )*dΩ
a_ref(dΩ) = a([f₁, f₂],uh,dΩ)
@test sum(a(uH,uH,dΩh)) ≈ sum(a(uh,uh,dΩh))
@test sum(a(uh,uH,dΩH)) ≈ sum(a(uH,uH,dΩH))
@test sum(a(uH,uH,dΩHh)) ≈ sum(a(uh,uh,dΩHh))

end
