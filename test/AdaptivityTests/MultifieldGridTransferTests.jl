module GridTransferTests

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
model = refine(cart_model; num_refinements=2)
trian = Triangulation(model)

# Triangulations
ftrian = Triangulation(get_model(model))
ctrian = Triangulation(get_parent(model))

# Measures
dΩ_f   = Measure(ftrian,2)
dΩ_c   = Measure(ctrian,2)
dΩ_cf  = Measure(ctrian,trian,1)

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

# Test values
pts = [VectorValue(rand(2)) for i=1:10]
uHi₁,uHi₂ = uHi
uhi₁,uhi₂ = uhi
for pt in pts
  @test uhi₁(pt) ≈ f₁(pt)
  @test uhi₂(pt) ≈ f₂(pt)
  @test uHi₁(pt) ≈ f₁(pt)
  @test uHi₂(pt) ≈ f₂(pt)
end

end
