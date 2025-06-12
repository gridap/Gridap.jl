module PolytopalHDGConvgTests

using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.Arrays

using Gridap.FESpaces: get_facet_diameter

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  P = LocalOperator(
    LocalSolveMap(), V0, mass, mass; trian_out = Ω
  )
  return P
end

##############################################################
#   HDG+ variant
##############################################################
function solve_Poisson_HDG(domain,nc,order,uex,f)
  model = UnstructuredDiscreteModel(simplexify(CartesianDiscreteModel(domain,nc)))
  vmodel = Gridap.Geometry.voronoi(model)

  D = num_cell_dims(vmodel)
  Ω = Triangulation(ReferenceFE{D}, vmodel)
  Γ = Triangulation(ReferenceFE{D-1}, vmodel)

  ptopo = Geometry.PatchTopology(vmodel)
  Ωp = Geometry.PatchTriangulation(vmodel,ptopo)
  Γp = Geometry.PatchBoundaryTriangulation(vmodel,ptopo)

  # HDG bulk spaces
  V = FESpaces.PolytopalFESpace(Ω, VectorValue{D, Float64}, order; space=:P)
  Q = FESpaces.PolytopalFESpace(Ω, Float64, order+1; space=:P)
  
  # HDG skeleton space(s)
  reffeM = ReferenceFE(lagrangian, Float64, order; space=:P) 
  M_test = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")
  M = TrialFESpace(M_test, uex)

  mfs = MultiField.BlockMultiFieldStyle(2,(2,1))
  X_full = MultiFieldFESpace([V, Q, M];style=mfs)
  X_elim = MultiFieldFESpace([V, Q])
  X_ret  = M

  polys = get_polytopes(vmodel)
  hT = map(p->FESpaces.get_facet_diameter(p,D),polys) 
  τT =  CellField(1 ./ hT,Ωp) # HDG stab parameter

  degree = 2*(order+1)
  dΩp = Measure(Ωp,degree)
  dΓp = Measure(Γp,degree)
  
  n   = get_normal_vector(Γp)
  Πn(u) = u⋅n
  Π(u) = change_domain(u,Γp,DomainStyle(u))
  
  PΓ = projection_operator(M_test, Γp, dΓp)   # Input has to be an UnconstrainedFESpace
  
  Pqhn(qh,uh,sh) = Πn(qh) + τT * (PΓ(uh)-sh) # Transmission condition
  a((qh,uh,sh),(vh,wh,lh)) = ∫( qh⋅vh - uh*(∇⋅vh) - qh⋅∇(wh) )dΩp + ∫(sh*Πn(vh))dΓp +
                              ∫(Pqhn(qh,uh,sh)*(Π(wh) + lh))dΓp 
  l((vh,wh,hatmh)) = ∫( f*wh )*dΩp
  
  op = MultiField.StaticCondensationOperator(ptopo,X_full,X_elim,X_ret,a,l)
  sh = solve(op.sc_op)
  qh, uh = MultiField.backward_static_condensation(op,sh)

  eu = uex - uh
  l2u = sqrt( sum( ∫( eu ⋅ eu )dΩp ) ) 

  h = maximum(hT)

  return l2u, h
end

function convg_test(domain,ncs,order,u,f)
  
  el2 = Float64[]
  hs  = Float64[]
  for nc in ncs
    l2u, h = solve_Poisson_HDG(domain,nc,order,u,f)
    println(l2u)
    push!(el2,l2u)
    push!(hs,h)
  end
  println(el2)
  println(hs)
  el2, hs
end

function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end


u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
q(x) = -∇(u)(x)
f(x) = (∇ ⋅ q)(x)

domain = (0,1,0,1)
ncs = [(2,2),(4,4),(8,8),(16,16),(32,32),(64,64),(128,128)]
order = 0

el, hs = convg_test(domain,ncs,order,u,f)
println("Slope L2-norm u: $(slope(hs,el))")

end # module