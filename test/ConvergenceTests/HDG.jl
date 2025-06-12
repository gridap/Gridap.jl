module HDGConvgTests

using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers

function solve_Poisson_HDG(domain,nc,order,uex,f)
  model = UnstructuredDiscreteModel(simplexify(CartesianDiscreteModel(domain,nc)))
  D = num_cell_dims(model)
  Ω = Triangulation(ReferenceFE{D}, model)
  Γ = Triangulation(ReferenceFE{D-1}, model)

  ptopo = Geometry.PatchTopology(model)
  Ωp = Geometry.PatchTriangulation(model,ptopo)
  Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

  # Reference FEs
  reffeV = ReferenceFE(lagrangian, VectorValue{D, Float64}, order; space=:P)
  reffeQ = ReferenceFE(lagrangian, Float64, order; space=:P)
  reffeM = ReferenceFE(lagrangian, Float64, order; space=:P)

  # HDG test FE Spaces
  V_test = TestFESpace(Ω, reffeV; conformity=:L2)
  Q_test = TestFESpace(Ω, reffeQ; conformity=:L2)
  M_test = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")

  # HDG trial FE Spaces
  V = TrialFESpace(V_test)
  Q = TrialFESpace(Q_test)
  M = TrialFESpace(M_test, uex)

  mfs = MultiField.BlockMultiFieldStyle(2,(2,1))
  X_full = MultiFieldFESpace([V, Q, M];style=mfs)
  X_elim = MultiFieldFESpace([V, Q])
  X_ret = M

  τ = 1.0 # HDG stab parameter

  degree = 2*(order+1)
  dΩp = Measure(Ωp,degree)
  dΓp = Measure(Γp,degree)
  n = get_normal_vector(Γp)

  Πn(u) = u⋅n
  Π(u) = change_domain(u,Γp,DomainStyle(u))
  a((qh,uh,sh),(vh,wh,lh)) = ∫( qh⋅vh - uh*(∇⋅vh) - qh⋅∇(wh) )dΩp + ∫(sh*Πn(vh))dΓp +
                            ∫((Πn(qh) + τ*(Π(uh) - sh))*(Π(wh) + lh))dΓp
  l((vh,wh,lh)) = ∫( f*wh )*dΩp

  op = MultiField.StaticCondensationOperator(ptopo,X_full,X_elim,X_ret,a,l)

  sh = solve(op.sc_op)
  qh, uh = MultiField.backward_static_condensation(op,sh)

  eu = u - uh
  l2u = sqrt( sum( ∫( eu ⋅ eu )dΩp ) ) 

  return l2u
end

function convg_test(domain,ncs,order,u,f)
  el2 = Float64[]
  hs = Float64[]
  for nc in ncs
    l2 = solve_Poisson_HDG(domain,nc,order,u,f)
    println(l2)
    push!(el2,l2)
    h = 1/nc[1]
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
ncs = [(2,2),(4,4),(8,8),(16,16),(32,32),(64,64),(128,128),(256,256)]
order = 1

el, hs = convg_test(domain,ncs,order,u,f)
println("Slope L2-norm u: $(slope(hs,el))")

end # module