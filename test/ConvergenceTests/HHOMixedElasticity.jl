module HHOMixxedElasticutyConvgTests

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs, Gridap.TensorValues
using Gridap.Arrays

ν(f) = skew_symmetric_gradient(f)

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  P = LocalOperator(
    LocalSolveMap(), V0, mass, mass; trian_out = Ω
  )
  return P
end

function gradient_operator(ptopo,X,order,Ω,Γp,dΩp,dΓp)
  reffe_L = ReferenceFE(lagrangian,SymTensorValue{2,Float64}, order; space = :P)
  L = FESpace(Ω, reffe_L; conformity=:L2)

  nΓ = get_normal_vector(Γp)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(q,))   =  ∫(p⊙q)dΩp
  rhs((uT,uF),(q,)) =  ∫(ε(uT)⊙q)dΩp + ∫((uF - Π(uT,Γp)) ⋅ (q⋅nΓ))dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  W = MultiFieldFESpace([L];style=MultiField.BlockMultiFieldStyle(1))
  G = LocalOperator(
    LocalSolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return G
end

function solve_Elasticity_HHO(domain,nc,order,uex,f,C)

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

  ##########################
  # Mixed order variant
  ##########################
  reffe_V = ReferenceFE(lagrangian, VectorValue{D,Float64}, order+1; space=:P)   # Bulk space
  reffe_M = ReferenceFE(lagrangian, VectorValue{D,Float64}, order; space=:P)     # Skeleton space
  V = FESpace(Ω, reffe_V; conformity=:L2)
  M = FESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")
  N = TrialFESpace(M,uex)

  mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
  X  = MultiFieldFESpace([V, N];style=mfs)
  Y  = MultiFieldFESpace([V, M];style=mfs)
  Xp = FESpaces.PatchFESpace(X,ptopo)

  PΓ = projection_operator(M, Γp, dΓp)
  G  = gradient_operator(ptopo,X,order,Ωp,Γp,dΩp,dΓp)

  global_assem = SparseMatrixAssembler(X,Y)
  patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

  function a(u,v)
    Gu_Ω, Gu_Γ = G(u)
    Gv_Ω, Gv_Γ = G(v)
    Gu = Gu_Ω + Gu_Γ
    Gv = Gv_Ω + Gv_Γ
    return ∫(Gu⊙C⊙Gv)dΩp
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

  function weakform()
    u, v = get_trial_fe_basis(X), get_fe_basis(Y)
    data1 = FESpaces.collect_cell_matrix_and_vector(X,Y,s(u,v),l(v),zero(X))
    data2 = FESpaces.collect_cell_matrix_and_vector(Xp,Xp,a(u,v),DomainContribution(),zero(Xp))
    data = FESpaces.merge_assembly_matvec_data(data1,data2)
    assemble_matrix_and_vector(global_assem,data)
  end

  function patch_weakform()
    u, v = get_trial_fe_basis(X), get_fe_basis(Y)
    data1 = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,X,Y,s(u,v),l(v),zero(X))
    data2 = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,Xp,a(u,v),DomainContribution(),zero(Xp))
    data = FESpaces.merge_assembly_matvec_data(data1,data2)
    return assemble_matrix_and_vector(patch_assem,data)
  end

  # Static condensation
  op = MultiField.StaticCondensationOperator(X,patch_assem,patch_weakform())
  ui, ub = solve(op)

  eu = uex - ui
  l2u = sqrt( sum(∫(eu ⋅ eu)dΩp) )
  h1u = l2u + sqrt( sum(∫(ε(eu) ⊙ ε(eu))dΩp) )

  h = maximum( sqrt(2)*sqrt.( get_array(∫(1)dΩp) ) )
  return l2u, h1u, h
end

function convg_test(domain,ncs,order,uex,f,C)
  el2 = Float64[]
  eh1 = Float64[]
  hs = Float64[]
  for nc in ncs
    println("  > Solving for nc = $(nc)")
    l2, h1, h = solve_Elasticity_HHO(domain,nc,order,uex,f,C)
    push!(el2,l2)
    push!(eh1,h1)
    push!(hs,h)
  end
  println("L2 error = ", el2)
  println("H1 error = ", eh1)
  el2, eh1, hs
end

function slope(hs,errors)
  x = log10.(hs)
  y = log10.(errors)
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end

μ = 1.0
λ = 1.0
I2 = one(SymTensorValue{2,Float64})
I4 = one(SymFourthOrderTensorValue{2,Float64})
C = λ * I2⊗I2 + μ * I4 # Isotropic elasticity tensor

u(x) = VectorValue(sin(2*π*x[1])*sin(2*π*x[2]),cos(2*π*x[1])*cos(2*π*x[2]))
g(x) = C⊙(ε(u)(x))
f(x) = -(∇⋅g)(x)

domain = (0,1,0,1)
ncs = [(2,2),(4,4),(8,8),(16,16),(32,32),(64,64),(128,128),(256,256)]
order = 1

el2s, eh1s, hs = convg_test(domain,ncs,order,u,f,C)
println("Slope L2-norm u: $(slope(hs,el2s))")
println("Slope H1-norm u: $(slope(hs,eh1s))")

end
