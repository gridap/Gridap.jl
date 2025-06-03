module HHOMixedConvgTests

using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs
using Gridap.Arrays

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  P = LocalOperator(
    LocalSolveMap(), V0, mass, mass; trian_out = Ω
  )
  return P
end

function reconstruction_operator(ptopo,L,X,Ω,Γp,dΩp,dΓp)
  reffe_Λ = ReferenceFE(lagrangian, Float64, 0; space=:P)
  Λ = FESpace(Ω, reffe_Λ; conformity=:L2)

  n = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅n
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ),(v,μ))   = ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
  rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = Y.multi_field_style
  W = MultiFieldFESpace([L,Λ];style=mfs)
  R = LocalOperator(
    LocalPenaltySolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end

##############################################################
#   HHO: Mixed order variant
##############################################################
function solve_Poisson_HHO(domain,nc,order,uex,f)
  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,nc))
  D = num_cell_dims(model)
  Ω = Triangulation(ReferenceFE{D}, model)
  Γ = Triangulation(ReferenceFE{D-1}, model)

  ptopo = Geometry.PatchTopology(model)
  Ωp = Geometry.PatchTriangulation(model,ptopo)
  Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

  qdegree = 2*(order+1)
  dΩp = Measure(Ωp,qdegree)
  dΓp = Measure(Γp,qdegree)

  reffe_V = ReferenceFE(lagrangian, Float64, order+1; space=:P)   # Bulk space
  reffe_M = ReferenceFE(lagrangian, Float64, order; space=:P)     # Skeleton space
  reffe_L = ReferenceFE(lagrangian, Float64, order+1; space=:P)   # Reconstruction space
  V = FESpace(Ω, reffe_V; conformity=:L2)
  M = FESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")
  L = FESpace(Ω, reffe_L; conformity=:L2)
  N = TrialFESpace(M,uex)

  mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
  X   = MultiFieldFESpace([V, N];style=mfs)
  Y   = MultiFieldFESpace([V, M];style=mfs) 
  Xp = FESpaces.PatchFESpace(X,ptopo)

  PΓ = projection_operator(M, Γp, dΓp)
  R  = reconstruction_operator(ptopo,L,Y,Ωp,Γp,dΩp,dΓp)

  patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

  function a(u,v)
    Ru_Ω, Ru_Γ = R(u)
    Rv_Ω, Rv_Γ = R(v)
    return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Ω) + ∇(Ru_Ω)⋅∇(Rv_Γ) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩp
  end

  hTinv =  CellField(1 ./ (sqrt(2).*sqrt.(get_array(∫(1)dΩp))),Ωp)
  function s(u,v)
    function SΓ(u)
      u_Ω, u_Γ = u
      return PΓ(u_Ω) - u_Γ
    end
    return ∫(hTinv * (SΓ(u)⋅SΓ(v)))dΓp
  end

  l((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp

  function patch_weakform()
    u, v = get_trial_fe_basis(X), get_fe_basis(Y)
    data1 = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,X,Y,s(u,v),l(v),zero(X))
    data2 = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,Xp,a(u,v),DomainContribution(),zero(Xp))
    data = FESpaces.merge_assembly_matvec_data(data1,data2)
    return assemble_matrix_and_vector(patch_assem,data)
  end

  # Static condensation
  op = MultiField.StaticCondensationOperator(X,V,N,patch_assem,patch_weakform())

  ub = solve(op.sc_op) 
  ui = MultiField.backward_static_condensation(op,ub)

  eu = uex - ui
  l2u= sqrt( sum(∫(eu ⋅ eu)dΩp) )
  h1u= l2u + sqrt( sum(∫(∇(eu) ⋅ ∇(eu))dΩp) )

  h = maximum( sqrt(2)*sqrt.( get_array(∫(1)dΩp) ) )

  return l2u, h1u, h
end

function convg_test(domain,ncs,order,uex,f)
    
  el2 = Float64[]
  eh1 = Float64[]
  hs = Float64[]
  for nc in ncs
    l2, h1, h = solve_Poisson_HHO(domain,nc,order,uex,f)
    push!(el2,l2)
    push!(eh1,h1)
    push!(hs,h)
  end
  println(el2)
  println(eh1)
  println(hs)
  el2, eh1, hs
end

function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end

uex(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
f(x) = -Δ(uex)(x)

domain = (0,1,0,1)
ncs = [(2,2),(4,4),(8,8),(16,16),(32,32),(64,64),(128,128),(256,256)]
order = 1

el2s, eh1s, hs = convg_test(domain,ncs,order,uex,f)
println("Slope L2-norm u: $(slope(hs,el2s))")
println("Slope H1-norm u: $(slope(hs,eh1s))")

end # module