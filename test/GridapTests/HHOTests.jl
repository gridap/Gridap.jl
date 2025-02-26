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
  
  W = MultiFieldFESpace([L,Λ];style=mfs)
  R = FESpaces.LocalOperator(
    FESpaces.HHO_ReconstructionOperatorMap(), ptopo, W, X, lhs, rhs; space_out = L
  )
  return R
end

function patch_dof_ids_w_bcs(ptopo,M::FESpace,MD::FESpace)
  ids_D = get_cell_dof_ids(MD)
  ids = get_cell_dof_ids(M)
  
  id_map = zeros(Int32,num_free_dofs(M))
  for (I,ID) in zip(ids,ids_D)
    id_map[I] = ID
  end
  
  patch_dofs = FESpaces.get_patch_dofs(M,ptopo)
  patch_dofs_D = lazy_map(Broadcasting(Reindex(id_map)),patch_dofs)

  free_values = zero_free_values(MD)
  dirichlet_values = get_dirichlet_dof_values(MD)
  cell_values = lazy_map(Broadcasting(Gridap.Arrays.PosNegReindex(free_values,dirichlet_values)),patch_dofs_D)
  cell_mask = collect(Bool,lazy_map(dofs -> any(dof -> dof < 0, dofs),patch_dofs_D))

  return patch_dofs_D, cell_values, cell_mask
end

function patch_dof_ids_w_bcs(ptopo,X::MultiFieldFESpace,XD::MultiFieldFESpace)
  nfields = MultiField.num_fields(X)

  pids, cv, cm = [], [], []
  map(X,XD) do M, MD
    a, b, c = patch_dof_ids_w_bcs(ptopo,M,MD)
    push!(pids,a)
    push!(cv,b)
    push!(cm,c)
  end

  patch_dofs = lazy_map(BlockMap(nfields,collect(1:nfields)),pids...)
  cell_values = lazy_map(BlockMap(nfields,collect(1:nfields)),cv...)
  cell_mask = map(any,zip(cm...))

  return patch_dofs, cell_values, cell_mask
end

##############################################################
u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
f(x) = -Δ(u)(x)

nc = (16,16)
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
M = FESpace(Γ, reffe_M; conformity=:L2)
L = FESpace(Ω, reffe_L; conformity=:L2)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
X   = MultiFieldFESpace([V, M];style=mfs)


PΓ = projection_operator(M, Γp, dΓp)
R = reconstruction_operator(ptopo,L,X,Ωp,Γp,dΩp,dΓp)

M0 = TestFESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")  # For assembly and imposing boundary conditions
MD = TrialFESpace(M0, u)

Y0 = MultiFieldFESpace([V, M0];style=mfs)
YD = MultiFieldFESpace([V, MD];style=mfs)

function a(u,v)
  Ru_Ω, Ru_Γ = R(u)
  Rv_Ω, Rv_Γ = R(v)
  return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Ω) + ∇(Ru_Ω)⋅∇(Rv_Γ) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩp
end

hTinv = 1 / CellField(get_array(∫(1)dΩp),Ωp)
function s(u,v)
  function S(u)
    u_Ω, u_Γ = u
    return PΓ(u_Ω) - u_Γ
  end
  return ∫(hTinv * (S(u)⋅S(v)))dΓp
end

l((vΩ,vΓ)) = ∫(f⋅vΩ)dΩp

function weakform(u,v)
  assem = SparseMatrixAssembler(YD,YD)

  data = FESpaces.collect_cell_matrix_and_vector(YD,YD,s(u,v),l(v),zero(YD))
  matvecdata, matdata, vecdata = data

  dofs_a, cellvals, cellmask = patch_dof_ids_w_bcs(ptopo,X,YD)
  cell_matvec = attach_dirichlet(get_array(a(u,v)),cellvals,cellmask)
  push!(matvecdata[1],cell_matvec)
  push!(matvecdata[2],dofs_a)
  push!(matvecdata[3],dofs_a)

  data = (matvecdata,matdata,vecdata)
  assemble_matrix_and_vector(assem,data)
end

function patch_weakform(u,v)
  assem = FESpaces.PatchAssembler(ptopo,YD,YD)

  data = FESpaces.collect_patch_cell_matrix_and_vector(assem,YD,YD,s(u,v),l(v),zero(YD))
  matvecdata, matdata, vecdata = data

  dofs_a, cellvals, cellmask = patch_dof_ids_w_bcs(ptopo,X,YD)
  cell_matvec = attach_dirichlet(get_array(a(u,v)),cellvals,cellmask)
  p = Geometry.get_patch_faces(ptopo,num_cell_dims(model))
  q = Base.OneTo(num_cells(model))

  rows_a = FESpaces.attach_patch_map(FESpaces.AssemblyStrategyMap{:rows}(assem.strategy),[dofs_a],[q])[1]
  cols_a = FESpaces.attach_patch_map(FESpaces.AssemblyStrategyMap{:cols}(assem.strategy),[dofs_a],[q])[1]
  push!(matvecdata[1],cell_matvec)
  push!(matvecdata[2],rows_a)
  push!(matvecdata[3],cols_a)
  push!(matvecdata[4],p)

  data = (matvecdata,matdata,vecdata)
  return assemble_matrix_and_vector(assem,data)
end

function sc_assembly(u,v)
  full_matvecs = patch_weakform(u,v)
  sc_matvecs = lazy_map(FESpaces.StaticCondensationMap(),full_matvecs)

  assem = FESpaces.PatchAssembler(ptopo,YD,YD)
  patch_rows = assem.strategy.array.array[2,2].patch_rows
  patch_cols = assem.strategy.array.array[2,2].patch_cols
  matvecdata = ([sc_matvecs,],[patch_rows,],[patch_cols,])
  matdata = ([],[],[]) # dummy matdata
  vecdata = ([],[],[]) # dummy vecdata
  data = (matvecdata, matdata, vecdata)
  A, b = assemble_matrix_and_vector(SparseMatrixAssembler(MD,M0),data)
  return A, b
end

function backward_sc(xb)
  assem = FESpaces.PatchAssembler(ptopo,YD,YD)
  patch_rows_bb = assem.strategy.array.array[2,2].patch_rows

  patchwise_xb = map(patch_rows_bb) do patch_row
    patchwise_xb = xb[patch_row]
    return patchwise_xb
  end

  full_matvecs = patch_weakform(get_trial_fe_basis(YD),get_fe_basis(YD))
  patchwise_xi_vec = lazy_map(Gridap.FESpaces.BackwardStaticCondensationMap(),full_matvecs,patchwise_xb)

  patch_rows_ii = assem.strategy.array.array[1,1].patch_rows
  patch_cols_ii = assem.strategy.array.array[1,1].patch_cols
  patchwise_xi_vecdata = ([patchwise_xi_vec,],[patch_rows_ii,],[patch_cols_ii,])

  xi = assemble_vector(SparseMatrixAssembler(V,V), patchwise_xi_vecdata)
  wh = FEFunction(V,xi)
  return wh
end

A, b = weakform(get_trial_fe_basis(YD),get_fe_basis(YD))
x = A \ b

Asc, bsc = sc_assembly(get_trial_fe_basis(YD),get_fe_basis(YD))
xb = Asc \ bsc

wh = backward_sc(xb)
xi = get_free_dof_values(wh)
