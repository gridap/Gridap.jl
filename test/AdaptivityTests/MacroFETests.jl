using Test
using Gridap
using Gridap.Adaptivity, Gridap.Geometry, Gridap.ReferenceFEs, Gridap.Arrays, Gridap.Helpers
using Gridap.CellData, Gridap.Fields, Gridap.FESpaces
using FillArrays

using Gridap.Fields

using Gridap.Adaptivity: num_subcells

function test_macro_reffe(model,fmodel,rrule,order)
  poly = get_polytope(rrule)
  reffe = LagrangianRefFE(Float64,poly,order)
  sub_reffes = Fill(reffe,num_subcells(rrule))
  macro_reffe = Adaptivity.MacroReferenceFE(rrule,sub_reffes)

  macro_quad  = Quadrature(poly,Adaptivity.CompositeQuadrature(),rrule,2*order)
  macro_quad_bis = Adaptivity.CompositeQuadrature(Quadrature(poly,4*order),rrule)
  ReferenceFEs.test_quadrature(macro_quad)
  ReferenceFEs.test_quadrature(macro_quad_bis)

  Ω = Triangulation(model)
  Ωf = Triangulation(fmodel)

  dΩ = Measure(Ω,2*order)
  dΩm = Measure(Ω,macro_quad)
  dΩfc = Measure(Ω,Ωf,2*order)

  V = FESpace(model,reffe)
  Vm = FESpace(model,macro_reffe)
  Vf  = FESpace(fmodel,reffe)
  @assert num_free_dofs(Vm) == num_free_dofs(Vf)

  # um and uf should be the same, but we cannot compare them directly. We 
  # project them into a common coarse space.
  u_exact(x) = cos(2π*x[1])*sin(2π*x[2])
  um = interpolate(u_exact,Vm)
  uf = interpolate(u_exact,Vf)

  M = assemble_matrix((u,v) -> ∫(u*v)dΩ,V,V)
  bf = assemble_vector(v -> ∫(v*uf)dΩfc,V)
  bm = assemble_vector(v -> ∫(v*um)dΩm,V)
  @test norm(bf-bm) < 1.e-10
  xf = M\bf
  xm = M\bm
  @test norm(xf-xm) < 1.e-10
end

order = 3

model1 = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(4,4)))
fmodel1 = refine(model1)
rrule1 = Adaptivity.RedRefinementRule(QUAD)
test_macro_reffe(model1,fmodel1,rrule1,order)

model2 = simplexify(model1)
fmodel2 = refine(model2,refinement_method="barycentric")
rrule2 = Adaptivity.BarycentricRefinementRule(TRI)
test_macro_reffe(model2,fmodel2,rrule2,order)
