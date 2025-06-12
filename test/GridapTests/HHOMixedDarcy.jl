
"""
References: 

- The hybrid high-order method for polytopal meshes. Di Pietro, Droniou. 2020
 
"""

using Test
using LinearAlgebra
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Gridap.ReferenceFEs
using Gridap.Arrays

function l2_error(uh,u,dΩ)
  eh = uh - u
  return sqrt(sum(∫(eh⋅eh)*dΩ))
end

function swap_field_ids(u, ids, nfields)
  @assert length(ids) == length(u)
  fields = map((ui,id) -> swap_field_ids(ui,id,nfields), u, ids)
  return MultiField.MultiFieldCellField(fields)
end

function swap_field_ids(u::MultiField.MultiFieldFEBasisComponent, id, nfields)
  return MultiField.MultiFieldFEBasisComponent(u.single_field, id, nfields)
end

function projection_operator(V, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  V0 = FESpaces.FESpaceWithoutBCs(V)
  P = LocalOperator(
    LocalSolveMap(), V0, mass, mass; trian_out = Ω
  )
  return P
end

function divergence_operator(ptopo,L,X,Ω,Γp,dΩp,dΓp)
  nrel = get_normal_vector(Γp)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(q,))   =  ∫(p⋅q)dΩp
  rhs((uT,uF),(q,)) =  ∫((∇⋅uT)⋅q)dΩp + ∫((uF - Π(uT,Γp))⋅nrel * q)dΓp
  
  Y = FESpaces.FESpaceWithoutBCs(X)
  mfs = MultiField.BlockMultiFieldStyle(1)
  W = MultiFieldFESpace([L];style=mfs)
  D = LocalOperator(
    LocalSolveMap(), ptopo, W, Y, lhs, rhs; space_out = L
  )
  _D(u) = swap_field_ids(D(swap_field_ids(u,[1,2],2)),[1,2],3)
  return _D
end

function hdiv_reconstruction_operator(ptopo, order, X, Ω, Γ, Γp, dΩp, dΓp)
  L = FESpaces.PolytopalFESpace(Ω, Float64, order; space=:RT)
  LT = FESpaces.PolytopalFESpace(Ω, VectorValue{2,Float64}, order-1; space=:P)
  LF = FESpaces.PolytopalFESpace(Γ, Float64, order; space=:P)

  nrel = get_normal_vector(Γp)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs((p,),(qT,qF)) = ∫(p⋅qT)dΩp + ∫(Π(p,Γp)⋅nrel * qF)dΓp
  rhs((uT,uF),(qT,qF)) = ∫(uT⋅qT)dΩp + ∫(uF⋅nrel * qF)dΓp

  W1 = MultiFieldFESpace([L]; style=MultiField.BlockMultiFieldStyle(1))
  W2 = MultiFieldFESpace([LT,LF]; style=MultiField.BlockMultiFieldStyle(1,(2,)))
  Y = FESpaces.FESpaceWithoutBCs(X)

  RT = LocalOperator(
    LocalSolveMap(RowMaximum()), ptopo, W1, Y, lhs, rhs; space_out = L, space_test = W2
  )
  _RT(u) = swap_field_ids(RT(swap_field_ids(u,[1,2],2)),[1,2],3)
  return _RT
end

##############################################################
# u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
u(x) = VectorValue(x[1],-x[2])
p(x) = x[1] + x[2]
f(x) = u(x) - ∇(p)(x)

# NOTE: Only works for simpolicial meshes, otherwise the RT-reconstruction is not well defined.
n = 8
model = simplexify(CartesianDiscreteModel((0,1,0,1),(n,n)))
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

order = 1
qdegree = 2*(order+1)

dΩ = Measure(Ω,qdegree)
dΩp = Measure(Ωp,qdegree)
dΓp = Measure(Γp,qdegree)

##########################
# Mixed order variant
##########################

reffe_V = ReferenceFE(lagrangian, VectorValue{2,Float64}, order+1; space=:P)   # Bulk space
reffe_M = ReferenceFE(lagrangian, VectorValue{2,Float64}, order; space=:P)     # Skeleton space
reffe_Q = ReferenceFE(lagrangian, Float64, order; space=:P) # Pressure space
V = FESpace(Ω, reffe_V; conformity=:L2)
M = FESpace(Γ, reffe_M; conformity=:L2, dirichlet_tags="boundary")
Q = FESpace(Ω, reffe_Q; conformity=:L2) # Pressure space
N = TrialFESpace(M,u)

mfs = MultiField.BlockMultiFieldStyle(2,(1,2))
X  = MultiFieldFESpace([V, N, Q];style=mfs)
Y  = MultiFieldFESpace([V, M, Q];style=mfs)
Xp = FESpaces.PatchFESpace(X,ptopo)

YR = MultiFieldFESpace([V, M];style=MultiField.BlockMultiFieldStyle(2,(1,1)))

PΓ = projection_operator(M, Γp, dΓp)
D  = divergence_operator(ptopo,Q,YR,Ωp,Γp,dΩp,dΓp)
R = hdiv_reconstruction_operator(ptopo, order, YR, Ω, Γ, Γp, dΩp, dΓp)

global_assem = SparseMatrixAssembler(X,Y)
patch_assem = FESpaces.PatchAssembler(ptopo,X,Y)

# Check that Div ∘ R = D
v = get_fe_basis(YR);
Rv = R(v) 
Dv = D(v)
vals = get_array(∫((∇⋅Rv[1] - Dv[1]) + (∇⋅Rv[2] - Dv[2]))dΩp)
map(vals) do vec
  all(arr -> all(x -> abs(x) < 1.e-10, arr), vec.array[1:2])
end |> all

function a(Ru,Rv)
  u_Ω, u_Γ = Ru
  v_Ω, v_Γ = Rv
  return ∫(u_Ω⋅v_Ω + u_Γ⋅v_Ω + u_Ω⋅v_Γ + u_Γ⋅v_Γ)dΩp
end

function b(Ru,q)
  Ru_Ω, Ru_Γ = Ru
  return ∫((∇⋅Ru_Ω + ∇⋅Ru_Γ)⋅q)dΩp
end

hTinv =  CellField(1 ./ (sqrt(2).*sqrt.(get_array(∫(1)dΩp))),Ωp)
function s(u,v)
  function SΓ(u)
    u_Ω, u_Γ = u
    return PΓ(u_Ω) - u_Γ
  end
  return ∫(hTinv * (SΓ(u)⋅SΓ(v)))dΓp
end

function l(Rv)
  v_Ω, v_Γ = Rv
  return ∫(f⋅v_Ω)dΩp
end

function weakform(x,y)
  u, p = (x[1],x[2]), x[3]
  v, q = (y[1],y[2]), y[3]

  Ru, Rv = R(u), R(v)
  data1 = FESpaces.collect_cell_matrix_and_vector(X,Y,s(u,v),DomainContribution(),zero(X))
  data2 = FESpaces.collect_cell_matrix_and_vector(Xp,Xp,a(Ru,Rv),l(Rv),zero(Xp))
  data3 = FESpaces.collect_cell_matrix_and_vector(Xp,X,b(Ru,q),DomainContribution(),zero(Xp))
  data4 = FESpaces.collect_cell_matrix_and_vector(X,Xp,b(Rv,p),DomainContribution(),zero(X))
  data = FESpaces.merge_assembly_matvec_data(data1,data2,data3,data4)
  assemble_matrix_and_vector(global_assem,data)
end

function patch_weakform()
  u, v = get_trial_fe_basis(X), get_fe_basis(Y)
  data1 = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,X,Y,s(u,v),l(v),zero(X))
  data2 = FESpaces.collect_patch_cell_matrix_and_vector(patch_assem,Xp,Xp,a(u,v),DomainContribution(),zero(Xp))
  data = FESpaces.merge_assembly_matvec_data(data1,data2)
  return assemble_matrix_and_vector(patch_assem,data)
end

x, y = get_trial_fe_basis(X), get_fe_basis(Y)

# Monolithic solve
A, B = weakform(x,y)
xh = FEFunction(X, A \ B)
uhT, uhF, ph = xh[1], xh[2], xh[3]

l2_uT = l2_error(uhT,u,dΩp)
l2_uF = l2_error(uhF,u,dΓp)
l2_p = l2_error(ph,p,dΩp)

# Static condensation
op = MultiField.StaticCondensationOperator(X,V,N,patch_assem,patch_weakform())

ub = solve(op.sc_op) 
ui = MultiField.backward_static_condensation(op,ub)

l2_ui = sqrt(sum(∫((ui - u)⋅(ui - u))*dΩ))
@test l2_ui < 1e-10
