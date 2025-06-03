module HHOPolytopalConvgTests

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

function mf_projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(q,)) = ∫(p⋅q)dΩ
  rhs((uT,uF),(q,)) = ∫(Π(uT,Ω)⋅q + Π(uF,Ω)⋅q)dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  W = MultiFieldFESpace([V0],style=MultiField.BlockMultiFieldStyle())
  P = LocalOperator(
    LocalSolveMap(), W, lhs, rhs; trian_out = Ω, space_out = V
  )
  return P
end

function patch_projection_operator(ptopo,V,X,Ω,dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(q,)) = ∫(p⋅q)dΩ
  rhs((uT,uF),(q,)) = ∫(Π(uT,Ω)⋅q + Π(uF,Ω)⋅q)dΩ

  V0 = FESpaces.FESpaceWithoutBCs(V)
  Y = FESpaces.FESpaceWithoutBCs(X)
  W = MultiFieldFESpace([V0],style=MultiField.BlockMultiFieldStyle())
  P = LocalOperator(
    LocalSolveMap(), ptopo, W, Y, lhs, rhs; trian_out = Ω, space_out = V
  )
  return P
end

function reconstruction_operator(ptopo,L,X,Ω,Γp,dΩp,dΓp)

  Λ = FESpaces.PolytopalFESpace(Ω, Float64, 0; space=:P)

  nrel = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅nrel
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ),(v,μ))   =  ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
  rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = Y.multi_field_style
  W = MultiFieldFESpace([L,Λ];style=mfs)
  R = LocalOperator(
    LocalPenaltySolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end

function potential_reconstruction(ptopo, X, L, R, uΩ, uΓ)
  u = get_trial_fe_basis(X);
  RuΩ, RuΓ =  R(u)
  
  Xp = FESpaces.PatchFESpace(X,ptopo)
  cvΩ = FESpaces.scatter_free_and_dirichlet_values(Xp[1],get_free_dof_values(uΩ),get_dirichlet_dof_values(X[1]))
  cvΓ = FESpaces.scatter_free_and_dirichlet_values(Xp[2],get_free_dof_values(uΓ),get_dirichlet_dof_values(X[2]))

  coeffs_Ω = map(*,get_data(RuΩ).args[1].args[1].args[1],cvΩ)
  coeffs_Γ = map(*,get_data(RuΓ).args[1].args[1].args[1],cvΓ)
  coeffs = map((a,b) -> a .+ b,coeffs_Ω,coeffs_Γ)

  return FEFunction(L,FESpaces.gather_free_values(L,coeffs))
end  

##############################################################

function solve_Poisson_HHO(domain,nc,order,uex,f)
  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,nc))
  vmodel = Gridap.Geometry.voronoi(simplexify(model))
  D = num_cell_dims(vmodel)
  Ω = Triangulation(ReferenceFE{D}, vmodel)
  Γ = Triangulation(ReferenceFE{D-1}, vmodel)

  ptopo = Geometry.PatchTopology(vmodel)
  Ωp = Geometry.PatchTriangulation(vmodel,ptopo)
  Γp = Geometry.PatchBoundaryTriangulation(vmodel,ptopo)

  qdegree = 2*(order+1)

  dΩp = Measure(Ωp,qdegree)
  dΓp = Measure(Γp,qdegree)

  reffe_M = ReferenceFE(lagrangian, Float64, order; space=:P)  # Skeleton space
  V = FESpaces.PolytopalFESpace(Ω, Float64, order; space=:P)   # Bulk Space
  M = FESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")
  L = FESpaces.PolytopalFESpace(Ω, Float64, order+1; space=:P) # Reconstruction space
  N = TrialFESpace(M,uex)

  mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
  X  = MultiFieldFESpace([V, N];style=mfs)
  Y  = MultiFieldFESpace([V, M];style=mfs)
  Xp = FESpaces.PatchFESpace(X,ptopo)

  R  = reconstruction_operator(ptopo,L,Y,Ωp,Γp,dΩp,dΓp)
  PΓ = projection_operator(M, Γp, dΓp)
  PΓ_mf = mf_projection_operator(M, Γp, dΓp)
  PΩ_mf = patch_projection_operator(ptopo,V,Xp,Ωp,dΩp)

  function a(Ru,Rv)
    Ru_Ω, Ru_Γ = Ru
    Rv_Ω, Rv_Γ = Rv
    return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Ω) + ∇(Ru_Ω)⋅∇(Rv_Γ) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩp
  end

  l((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp

  hFinv =  CellField(1 ./ get_array(∫(1)dΓp),Γp)
  function SΓa(u)
    u_Ω, u_Γ = u
    return PΓ(u_Ω) - u_Γ 
  end
  function SΓb(Ru)
    PΓRu_Ω, PΓRu_Γ = PΓ_mf(Ru)
    PΓPΩRu_Ω, PΓPΩRu_Γ = PΓ_mf(PΩ_mf(Ru))
    return (PΓRu_Ω - PΓPΩRu_Ω) + (PΓRu_Γ - PΓPΩRu_Γ)
  end

  patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

  function patch_weakform(u,v)
    mass_Γ(u,v) = ∫(hFinv*(u*v))dΓp
    Ru, Rv = R(u), R(v)
    SΓa_u, SΓa_v = SΓa(u), SΓa(v)
    SΓb_u, SΓb_v = SΓb(Ru), SΓb(Rv)
    aa = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,X,X,mass_Γ(SΓa_u,SΓa_v),l(v),zero(X))
    ab = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,X,Xp,mass_Γ(SΓa_u,SΓb_v),DomainContribution(),zero(X))
    ba = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,X,mass_Γ(SΓb_u,SΓa_v),DomainContribution(),zero(Xp))
    bb = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,Xp,mass_Γ(SΓb_u,SΓb_v),DomainContribution(),zero(Xp))
    ll = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,Xp,a(Ru,Rv),DomainContribution(),zero(Xp))
    data = FESpaces.merge_assembly_matvec_data(aa,ab,ba,bb,ll)
    return assemble_matrix_and_vector(patch_assem,data)
  end

  # Static condensation
  u, v = get_trial_fe_basis(X), get_fe_basis(Y);
  op = MultiField.StaticCondensationOperator(X,V,N,patch_assem,patch_weakform(u,v))

  uΓ = solve(op.sc_op) 
  uΩ = MultiField.backward_static_condensation(op,uΓ)

  cvΩ = FESpaces.scatter_free_and_dirichlet_values(Xp[1],get_free_dof_values(uΩ),get_dirichlet_dof_values(X[1]))
  cvΓ = FESpaces.scatter_free_and_dirichlet_values(Xp[2],get_free_dof_values(uΓ),get_dirichlet_dof_values(X[2]))

  RuΩ, RuΓ =  R(u)
  coeffs_Ω = map(*,get_data(RuΩ).args[1].args[1].args[1],cvΩ)
  coeffs_Γ = map(*,get_data(RuΓ).args[1].args[1].args[1],cvΓ)
  coeffs = map((a,b) -> a .+ b,coeffs_Ω,coeffs_Γ)

  Ruh = potential_reconstruction(ptopo, X, L, R, uΩ, uΓ)
  eRuh = Ruh - uex
  l2u = sqrt(sum( ∫(eRuh * eRuh)dΩp))
  h1u = l2u + sqrt(sum( ∫(∇(eRuh) ⋅ ∇(eRuh))dΩp))

  # eu  = uΩ - uex 
  # l2u = sqrt(sum( ∫(eu * eu)dΩp))
  # h1u = l2u + sqrt(sum( ∫(∇(eu) ⋅ ∇(eu))dΩp))

  polys = get_polytopes(vmodel)
  h = maximum( map(x -> FESpaces.get_facet_diameter(x,D), polys) )

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
order = 0

el2s, eh1s, hs = convg_test(domain,ncs,order,uex,f)
println("Slope L2-norm u: $(slope(hs,el2s))")
println("Slope H1-norm u: $(slope(hs,eh1s))")

end # module