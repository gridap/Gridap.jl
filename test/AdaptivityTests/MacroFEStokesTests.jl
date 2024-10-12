module MacroFEStokesTests
using Test
using Gridap
using Gridap.Adaptivity, Gridap.Geometry, Gridap.ReferenceFEs, Gridap.FESpaces, Gridap.Arrays
using Gridap.ReferenceFEs

using FillArrays

function main(Dc,reftype)
  @assert Dc ∈ [2,3]
  @assert reftype ∈ [:barycentric,:powellsabin]

  u_sol(x) = (Dc == 2) ? VectorValue(x[1],-x[2]) : VectorValue(x[1],-x[2],0.0)
  p_sol(x) = x[1] - 1.0/2.0

  domain = (Dc == 2) ? (0,1,0,1) : (0,1,0,1,0,1)
  nc = (Dc == 2) ? (2,2) : (1,1,1)
  model = simplexify(CartesianDiscreteModel(domain,nc))

  poly  = (Dc == 2) ? TRI : TET
  order = (reftype == :barycentric) ? Dc : Dc-1
  rrule = (reftype == :barycentric) ? Adaptivity.BarycentricRefinementRule(poly) : Adaptivity.PowellSabinRefinementRule(poly)

  subreffes_u = Fill(LagrangianRefFE(VectorValue{Dc,Float64},poly,order),Adaptivity.num_subcells(rrule))
  reffe_u = Adaptivity.MacroReferenceFE(rrule,subreffes_u)
  
  subreffes_p = Fill(LagrangianRefFE(Float64,poly,order-1),Adaptivity.num_subcells(rrule))
  reffe_p = Adaptivity.MacroReferenceFE(rrule,subreffes_p;conformity=L2Conformity())

  qdegree = 2*(order-1)
  quad  = Quadrature(poly,Adaptivity.CompositeQuadrature(),rrule,qdegree)

  V = FESpace(model,reffe_u,dirichlet_tags=["boundary"])
  Q = FESpace(model,reffe_p,conformity=:L2,constraint=:zeromean)
  U = TrialFESpace(V,u_sol)

  Ω = Triangulation(model)
  dΩ = Measure(Ω,quad)
  @assert abs(sum(∫(p_sol)dΩ)) < 1.e-15
  @assert abs(sum(∫(divergence(u_sol))dΩ)) < 1.e-15

  X = MultiFieldFESpace([U,Q])
  Y = MultiFieldFESpace([V,Q])

  f(x) = -Δ(u_sol)(x) + ∇(p_sol)(x)
  lap(u,v) = ∫(∇(u)⊙∇(v))dΩ
  a((u,p),(v,q)) = lap(u,v) + ∫(divergence(u)*q - divergence(v)*p)dΩ
  l((v,q)) = ∫(f⋅v)dΩ

  op = AffineFEOperator(a,l,X,Y)
  xh = solve(op)

  uh, ph = xh
  eh_u = uh - u_sol
  eh_p = ph - p_sol
  @test sum(∫(eh_u⋅eh_u)dΩ) < 1.e-10
  if reftype != :powellsabin
    @test sum(∫(eh_p*eh_p)dΩ) < 1.e-10
  end
end

# NOTE: Powell-Sabin split not working yet. The issue is we would need a global cellmap 
# directly from the sub-cells to the physical domain (due to how the split is built). 
# This is something I may do in the future. 

main(2,:barycentric)
#main(2,:powellsabin)
main(3,:barycentric)
#main(3,:powellsabin)

end # module