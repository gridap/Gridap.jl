using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  P = FESpaces.LocalOperator(
    FESpaces.LocalSolveMap(), V, mass, mass; trian_out = Ω
  )
  return P
end

function reconstruction_operator(ptopo,L,X,Ωp,Γp,dΩp,dΓp)
  reffe_Λ = ReferenceFE(lagrangian, Float64, 0; space=:P)
  Λ = FESpace(Ω, reffe_Λ; conformity=:L2)

  nrel = get_normal_vector(Γp)
  Πn(v) = ∇(v)⋅nrel
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((u,λ),(v,μ))   = ∫( (∇(u)⋅∇(v)) + (μ*u) + (λ*v) )dΩp
  rhs((uT,uF),(v,μ)) =  ∫( (∇(uT)⋅∇(v)) + (uT*μ) )dΩp + ∫( (uF - Π(uT,Γp))*(Πn(v)) )dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  W = MultiFieldFESpace([L,Λ];style=mfs)
  R = FESpaces.LocalOperator(
    FESpaces.HHO_ReconstructionOperatorMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  return R
end

##############################################################
u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
f(x) = -Δ(u)(x)

nc = (2,2)
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),nc))
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

order = 1
qdegree = 2*(order+1)

dΩp = Measure(Ωp,qdegree)
dΓp = Measure(Γp,qdegree)

reffe_V = ReferenceFE(lagrangian, Float64, order+1; space=:P)   # Bulk space
reffe_M = ReferenceFE(lagrangian, Float64, order; space=:P)     # Skeleton space
reffe_L = ReferenceFE(lagrangian, Float64, order+1; space=:P)   # Reconstruction space
V = FESpace(Ω, reffe_V; conformity=:L2)
M = FESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")
L = FESpace(Ω, reffe_L; conformity=:L2)
N = TrialFESpace(M,u)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
X  = MultiFieldFESpace([V, N];style=mfs)
Y  = MultiFieldFESpace([V, M];style=mfs)
Xp = FESpaces.PatchFESpace(X,ptopo)

PΓ = projection_operator(M, Γp, dΓp)
R = reconstruction_operator(ptopo,L,Y,Ωp,Γp,dΩp,dΓp)

function a(u,v)
  Ru_Ω, Ru_Γ = R(u)
  Rv_Ω, Rv_Γ = R(v)
  return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Ω) + ∇(Ru_Ω)⋅∇(Rv_Γ) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩp
end

hTinv = 1 ./ CellField((sqrt(2).*sqrt.(get_array(∫(1)dΩp))),Ωp)
function s(u,v)
  function S(u)
    u_Ω, u_Γ = u
    return PΓ(u_Ω) - u_Γ
  end
  return ∫(hTinv * (S(u)⋅S(v)))dΓp
end

l((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp

function weakform(u,v)
  assem = SparseMatrixAssembler(X,Y)

  data1 = FESpaces.collect_cell_matrix_and_vector(X,Y,s(u,v),l(v),zero(X))
  data2 = FESpaces.collect_cell_matrix_and_vector(Xp,Xp,a(u,v),DomainContribution(),zero(Xp))
  data = FESpaces.merge_assembly_matvec_data(data1,data2)

  assemble_matrix_and_vector(assem,data)
end

function patch_weakform(u,v)
  assem = FESpaces.PatchAssembler(ptopo,X,Y)

  data1 = FESpaces.collect_patch_cell_matrix_and_vector(assem,X,Y,s(u,v),l(v),zero(X))
  data2 = FESpaces.collect_patch_cell_matrix_and_vector(assem,Xp,Xp,a(u,v),DomainContribution(),zero(Xp))
  data = FESpaces.merge_assembly_matvec_data(data1,data2)

  return assemble_matrix_and_vector(assem,data)
end

function sc_assembly(u,v)
  full_matvecs = patch_weakform(u,v)
  sc_matvecs = lazy_map(FESpaces.StaticCondensationMap(),full_matvecs)

  assem = FESpaces.PatchAssembler(ptopo,X,Y)
  patch_rows = assem.strategy.array.array[2,2].patch_rows
  patch_cols = assem.strategy.array.array[2,2].patch_cols
  matvecdata = ([sc_matvecs,],[patch_rows,],[patch_cols,])
  matdata = ([],[],[]) # dummy matdata
  vecdata = ([],[],[]) # dummy vecdata
  data = (matvecdata, matdata, vecdata)
  A, b = assemble_matrix_and_vector(SparseMatrixAssembler(N,M),data)
  return A, b
end

function backward_sc(xb)
  assem = FESpaces.PatchAssembler(ptopo,X,Y)
  patch_rows_bb = assem.strategy.array.array[2,2].patch_rows

  patchwise_xb = map(patch_rows_bb) do patch_row
    patchwise_xb = xb[patch_row]
    return patchwise_xb
  end

  full_matvecs = patch_weakform(get_trial_fe_basis(X),get_fe_basis(Y))
  patchwise_xi_vec = lazy_map(Gridap.FESpaces.BackwardStaticCondensationMap(),full_matvecs,patchwise_xb)

  patch_rows_ii = assem.strategy.array.array[1,1].patch_rows
  patch_cols_ii = assem.strategy.array.array[1,1].patch_cols
  patchwise_xi_vecdata = ([patchwise_xi_vec,],[patch_rows_ii,],[patch_cols_ii,])

  xi = assemble_vector(SparseMatrixAssembler(V,V), patchwise_xi_vecdata)
  wh = FEFunction(V,xi)
  return wh
end

A, b = weakform(get_trial_fe_basis(X),get_fe_basis(Y))
x = A \ b

Asc, bsc = sc_assembly(get_trial_fe_basis(X),get_fe_basis(Y))
xb = Asc \ bsc

wh = backward_sc(xb)
xi = get_free_dof_values(wh)
