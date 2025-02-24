using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs

function consistency(rec_vΩ, rec_uΩ, rec_vΓ, rec_uΓ, X)

  nfields = length(X.spaces)
  RvΩ = MultiField.MultiFieldFEBasisComponent(rec_vΩ,1,nfields)
  RuΩ = MultiField.MultiFieldFEBasisComponent(rec_uΩ,1,nfields)
  RvΓ = MultiField.MultiFieldFEBasisComponent(rec_vΓ,2,nfields)
  RuΓ = MultiField.MultiFieldFEBasisComponent(rec_uΓ,2,nfields)

    return ∫( (RuΩ*RvΩ)+(RuΩ*RvΓ)+(RuΓ*RvΩ)+(RuΓ*RvΓ) )Measure(Ω,4)
end

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

  return patch_dofs_D
end

function patch_dof_ids_w_bcs(ptopo,X::MultiFieldFESpace,XD::MultiFieldFESpace)
  nfields = MultiField.num_fields(X)
  pids = map(X,XD) do M, MD
    patch_dof_ids_w_bcs(ptopo,M,MD)
  end
  lazy_map(BlockMap(nfields,collect(1:nfields)),pids...)
end

##############################################################
u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
q(x) = -∇(u)(x)
f(x) = (∇ ⋅ q)(x)

nc = (2,2)
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),nc))
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

order = 0
qdegree = order+1

dΩ = Measure(Ω,qdegree)
dΩp = Measure(Ωp,qdegree)
dΓp = Measure(Γp,qdegree)

reffe_V = ReferenceFE(lagrangian, Float64, order; space=:P)   # Bulk space
reffe_M = ReferenceFE(lagrangian, Float64, order; space=:P)   # Skeleton space
reffe_L = ReferenceFE(lagrangian, Float64, order+1; space=:P) # Reconstruction space
V = FESpace(Ω, reffe_V; conformity=:L2)
M = FESpace(Γ, reffe_M; conformity=:L2)
L = FESpace(Ω, reffe_L; conformity=:L2)

mfs = MultiField.BlockMultiFieldStyle(2,(1,1))
X   = MultiFieldFESpace([V, M];style=mfs)

PΩ = projection_operator(V, Ω, dΩ)
PΓ = projection_operator(M, Γp, dΓp)

R = reconstruction_operator(ptopo,L,X,Ωp,Γp,dΩp,dΓp)

M0 = TestFESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")  # For assembly and imposing boundary conditions
MD = TrialFESpace(M0, u)

Y0 = MultiFieldFESpace([V, M0];style=mfs)
YD = MultiFieldFESpace([V, MD];style=mfs)

function a(u,v)
  Ru_Ω, Ru_Γ = R(u)
  Rv_Ω, Rv_Γ = R(v)
  return ∫(∇(Ru_Ω)⋅∇(Rv_Ω) + ∇(Ru_Γ)⋅∇(Rv_Γ))dΩ
end

hF = 1 / CellField(get_array(∫(1)dΓp),Γp)
function s(u,v)
  function S(u)
    u_Ω, u_Γ = u
    return PΓ(u_Ω) - u_Γ
  end
  return ∫(hF * (S(u)⋅S(v)))dΓp
end

function biform(u,v)
  c_a = a(u,v)
  c_s = s(u,v)

  dofs_a = patch_dof_ids_w_bcs(ptopo,X,YD)
  dofs_s = get_cell_dof_ids(YD,Γp)
  matdata = ([get_array(c_a),get_array(c_s)],[dofs_a,dofs_s],[dofs_a,dofs_s])
  assemble_matrix(SparseMatrixAssembler(YD,YD),matdata)
end

x = get_trial_fe_basis(YD)
y = get_fe_basis(YD)

A = biform(x,y)

