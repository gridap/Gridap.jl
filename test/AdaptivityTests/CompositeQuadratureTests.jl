module CompositeQuadratureTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

sol(x) = x[1] + x[2]

model1 = CartesianDiscreteModel((0,1,0,1),(4,4))
model2 = simplexify(model1)
model3 = CartesianDiscreteModel((0,1,0,1,0,1),(4,4,4))
model4 = simplexify(model3)
models = [model1,model2,model3,model4]

for (k,poly) in enumerate([QUAD,TRI,HEX,TET])
  model = models[k]
  ncells = num_cells(model)

  order = 1
  reffe = LagrangianRefFE(Float64,poly,order)
  V     = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
  U     = TrialFESpace(V,sol)

  rr1    = RefinementRule(poly,2)
  rr2    = RefinementRule(reffe,3)
  rrules = CompressedArray([rr1,rr2],rand([1,2],ncells))

  quad_order = 2*order + 1
  Ω = Triangulation(model)
  dΩ = Measure(Ω,quad_order)
  dΩ_rr = Measure(Ω,rrules,quad_order)

  u = interpolate(sol,U)
  a(u,v) = ∫(v⋅u)*dΩ
  a_rr(u,v) = ∫(v⋅u)*dΩ_rr

  @test sum(a_rr(1,1)) ≈ 1.0
  @test sum(a_rr(u,u)) ≈ sum(a(u,u))
end

end