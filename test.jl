using Test

using Gridap
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Geometry
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.MultiField

using ForwardDiff
using SparseArrays

domain = (0,1,0,1)
partition = (5,5)
model = CartesianDiscreteModel(domain,partition)
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model,tags=["tag_5"])
dΩ = Measure(Ω,2)
dΓ = Measure(Γ,2)
V1 = FESpace(Γ,ReferenceFE(lagrangian,Float64,1))
V2 = FESpace(model,ReferenceFE(lagrangian,VectorValue{2,Float64},1))
V3 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
X = MultiFieldFESpace([V1,V2,V3])
f(xh) = ∫(xh[1]+xh[2]⋅xh[2]+xh[1]*xh[3])dΓ
uh = zero(X)
du = gradient(f,uh)
du1 = gradient(x->f((x,uh[2],uh[3])),uh[1])
du2 = gradient(x->f((uh[1],x,uh[3])),uh[2])
du3 = gradient(x->f((uh[1],uh[2],x)),uh[3])

@test lazy_map(Gridap.MultiField.GetIndex(1),du[Γ]) == du1[Γ]
@test lazy_map(Gridap.MultiField.GetIndex(2),du[Γ]) == du2[Γ]
@test lazy_map(Gridap.MultiField.GetIndex(3),du[Γ]) == du3[Γ]

du1_vec = assemble_vector(du1,V1)
du2_vec = assemble_vector(du2,V2)
du3_vec = assemble_vector(du3,V3)
du_vec = assemble_vector(du,X)

@test du_vec == [du1_vec;du2_vec;du3_vec]

f2(xh,yh) = ∫(xh[1]⋅yh[1]+xh[2]⋅yh[2]+xh[1]⋅xh[2]⋅yh[2]+xh[1]*xh[3]*yh[3])dΓ
dv = get_fe_basis(X)
j = jacobian(uh->f2(uh,dv),uh)
J = assemble_matrix(j,X,X)

f2_jac(xh,dxh,yh) = ∫(dxh[1]⋅yh[1]+dxh[2]⋅yh[2]+dxh[1]⋅xh[2]⋅yh[2]+xh[1]⋅dxh[2]⋅yh[2]+dxh[1]*xh[3]*yh[3]+xh[1]*dxh[3]*yh[3])dΓ
op = FEOperator(f2,f2_jac,X,X)
J_fwd = jacobian(op,uh)

@test J_fwd == J

Λ = SkeletonTriangulation(model)
dΛ = Measure(Λ,2)
f(xh) = ∫(mean(xh[1])+mean(xh[2])⋅mean(xh[2])+mean(xh[1])*mean(xh[3]))dΛ
uh = zero(X)
du = gradient(f,uh)
du1 = gradient(x->f((x,uh[2],uh[3])),uh[1])
du2 = gradient(x->f((uh[1],x,uh[3])),uh[2])
du3 = gradient(x->f((uh[1],uh[2],x)),uh[3])

Xbis = MultiFieldFESpace([V1,V1,V1])
uhbis = zero(Xbis)
du_bis = gradient(f,uhbis;ad_type=:monolithic)

@test lazy_map(Gridap.MultiField.GetIndex(1),du[Λ]) == du1[Λ]
@test lazy_map(Gridap.MultiField.GetIndex(2),du[Λ]) == du2[Λ]
@test lazy_map(Gridap.MultiField.GetIndex(3),du[Λ]) == du3[Λ]

V4 = FESpace(Λ,ReferenceFE(lagrangian,Float64,1))

get_cell_dof_ids(V1,Λ)
get_cell_dof_ids(V2,Λ)
get_cell_dof_ids(V3,Λ)
get_cell_dof_ids(V4,Λ)

get_cell_dof_ids(X,Λ)

du1_vec = assemble_vector(du1,V1)
du2_vec = assemble_vector(du2,V2)
du3_vec = assemble_vector(du3,V3)
du_vec = assemble_vector(du,X)

@test du_vec == [du1_vec;du2_vec;du3_vec]

f2(xh,yh) = ∫(mean(xh[1])⋅mean(yh[1])+mean(xh[2])⋅mean(yh[2])+mean(xh[1])⋅mean(xh[2])⋅mean(yh[2])+mean(xh[1])*mean(xh[3])*mean(yh[3]))dΛ
dv = get_fe_basis(X)
j = jacobian(uh->f2(uh,dv),uh)
J = assemble_matrix(j,X,X)

f2_jac(xh,dxh,yh) = ∫(mean(dxh[1])⋅mean(yh[1])+mean(dxh[2])⋅mean(yh[2])+mean(dxh[1])⋅mean(xh[2])⋅mean(yh[2]) +
  mean(xh[1])⋅mean(dxh[2])⋅mean(yh[2])+mean(dxh[1])*mean(xh[3])*mean(yh[3])+mean(xh[1])*mean(dxh[3])*mean(yh[3]))dΛ
op = FEOperator(f2,f2_jac,X,X)
J_fwd = jacobian(op,uh)

@test J_fwd == J

# ## Hessian
# using ForwardDiff
# e((uh,ph,sh)) = ∫( uh*uh + uh*(ph ⋅ ph) + ph ⋅ ph)dΓ


# jacobian(k->gradient(s->e((s,uh[2],uh[3])),k),uh[1])


# # [d^2e/du^2 d^2e/dudp]
# # [d^2e/dpdu e^2e/dp^2 ]
# uh = interpolate([x->x[1],x->VectorValue(x[2],x[1]),x->x[2]],X)

# g_full = _change_argument(hessian,e,Γ,uh)
# g = map(i->_change_argument(hessian,Gridap.MultiField.restrict_function(e,uh,i),Γ,uh[i]),1:3)
# cell_u = get_cell_dof_values.(uh)
# cell_id = _compute_cell_ids.(uh,(Γ,))

# i_to_x = cell_u[1]
# j_to_i = cell_id[1]
# agrad = i_to_y -> Gridap.Arrays.autodiff_array_gradient(g[1],cell_u[1],cell_id[1])
# Gridap.Arrays.autodiff_array_jacobian.(agrad,cell_u[2],cell_id[2])

# # gradient
# # i_to_x = cell_u[1]
# # j_to_i = cell_id[1]

# dummy_tag = ()->()
# i_to_cfg_1 = lazy_map(FESpaces.ConfigMap(ForwardDiff.gradient,dummy_tag),cell_u[1])
# i_to_cfg_2 = lazy_map(FESpaces.ConfigMap(ForwardDiff.gradient,dummy_tag),cell_u[2])
# i_to_cfg_3 = lazy_map(FESpaces.ConfigMap(ForwardDiff.gradient,dummy_tag),cell_u[3])
# i_to_xdual_1 = lazy_map(FESpaces.DualizeMap(),i_to_cfg_1,cell_u[1])
# i_to_xdual_2 = lazy_map(FESpaces.DualizeMap(),i_to_cfg_2,cell_u[2])
# i_to_xdual_3 = lazy_map(FESpaces.DualizeMap(),i_to_cfg_3,cell_u[3])

# g_i_to_cfg2 = lazy_map(FESpaces.ConfigMap(ForwardDiff.gradient,dummy_tag),cell_u[2])
# g_i_to_xdual2 = lazy_map(FESpaces.DualizeMap(),g_i_to_cfg2,cell_u[2])

# # lazy_map(BlockMap(3,collect(Base.OneTo(3))),i_to_xdual_1,i_to_xdual_2,cell_u[3])
# # lazy_map((x...)->ArrayBlock([x[1],x[2],x[3]],[true,true,true]),i_to_xdual_1,i_to_xdual_2,cell_u[3])

# g_j_to_ydual = g_full(i_to_xdual_1,g_i_to_xdual2,cell_u[3])

# # g_j_to_ydual = g_full(i_to_xdual_1,i_to_xdual_2,cell_u[3])
# # g_j_to_ydual = g[2](g_i_to_xdual2)
# g_j_to_cfg = Gridap.Arrays.autodiff_array_reindex(i_to_cfg_1,cell_id[1])
# g_j_to_result = lazy_map(Gridap.Arrays.AutoDiffMap(),g_j_to_cfg,g_j_to_ydual)

# # # g_j_to_result = Gridap.Arrays.autodiff_array_gradient(g[2],cell_u[2],cell_id[2])

# # # _mf_cell_grad = lazy_map(BlockMap(3,collect(Base.OneTo(3))),i_to_xdual,cell_u[2],cell_u[3])
# # # _mf_cell_grad = lazy_map((x...)->ArrayBlock([x[1],x[2],x[3]],[true,true,true]),i_to_xdual,cell_u[2],cell_u[3])
# # # g_i_to_cfg = lazy_map(FESpaces.ConfigMap(ForwardDiff.gradient,dummy_tag),_mf_cell_grad)

# # g_i_to_cfg1 = lazy_map(FESpaces.ConfigMap(ForwardDiff.gradient,dummy_tag),i_to_xdual_1)
# # g_i_to_cfg2 = lazy_map(FESpaces.ConfigMap(ForwardDiff.gradient,dummy_tag),i_to_xdual_2)
# # g_i_to_cfg3 = lazy_map(FESpaces.ConfigMap(ForwardDiff.gradient,dummy_tag),i_to_xdual_3)

# # g_i_to_xdual1 = lazy_map(FESpaces.DualizeMap(),g_i_to_cfg1,i_to_xdual)
# # g_i_to_xdual2 = lazy_map(FESpaces.DualizeMap(),g_i_to_cfg2,cell_u[2])
# # g_i_to_xdual3 = lazy_map(FESpaces.DualizeMap(),g_i_to_cfg3,cell_u[3])

# # # full_g_i_to_xdual = lazy_map((x...)->ArrayBlock([x[1],x[2],x[3]],[true,true,true]),g_i_to_xdual2,g_i_to_xdual2,g_i_to_xdual3)

# # # g_j_to_ydual = g_full(_mf_cell_grad)

# # g_j_to_ydual = g[2](g_i_to_xdual2)

# # # g_j_to_ydual = g[1](g_i_to_xdual)
# # g_j_to_cfg = Gridap.Arrays.autodiff_array_reindex(g_i_to_cfg2,cell_id[2])
# # g_j_to_result = lazy_map(Gridap.Arrays.AutoDiffMap(),g_j_to_cfg,g_j_to_ydual)


# j_to_cfg = Gridap.Arrays.autodiff_array_reindex(i_to_cfg,j_to_i)
# j_to_result = lazy_map(Gridap.Arrays.AutoDiffMap(),j_to_cfg,g_j_to_result)
# j_to_result

# # agrad = i_to_y -> FESpaces.autodiff_array_gradient(g,i_to_y,cell_id[3])
# # cell_grad = FESpaces.autodiff_array_jacobian(agrad,cell_u[3],cell_id[3])

# # function autodiff_array_gradient_new(a,i_to_x,j_to_i)
# #   dummy_tag = ()->()
# #   i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient,dummy_tag),i_to_x)
# #   i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
# #   j_to_ydual = a(i_to_xdual)
# #   j_to_cfg = autodiff_array_reindex(i_to_cfg,j_to_i)
# #   j_to_result = lazy_map(AutoDiffMap(),j_to_cfg,j_to_ydual)
# #   j_to_result
# # end

# function autodiff_array_jacobian_new(a,i_to_x,j_to_i)
#   dummy_tag = ()->()
#   i_to_cfg = lazy_map(ConfigMap(ForwardDiff.jacobian,dummy_tag),i_to_x)
#   i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
#   j_to_ydual = a(i_to_xdual)
#   j_to_cfg = autodiff_array_reindex(i_to_cfg,j_to_i)
#   j_to_result = lazy_map(AutoDiffMap(),j_to_cfg,j_to_ydual)
#   j_to_result
# end

# # jacobian(x->gradient(e,x),uh)

# # function FESpaces._change_argument(op,f,trian,uh::GenericCellField)
# #   # U = get_fe_space(uh)
# #   function g(cell_u)
# #     cf = uh
# #     cell_grad = f(cf)
# #     get_contribution(cell_grad,trian)
# #   end
# #   g
# # end
# # FESpaces.get_cell_dof_values(uh::GenericCellField) = get_data(uh)
# # CellField(X,get_cell_dof_values(uh))
# # jacobian(x->gradient(e,x),uh)



# dv = get_fe_basis(X)

# terms = DomainContribution()
# F1 = Gridap.MultiField.restrict_function(e,uh,1)
# g1 = _change_argument(hessian,F1,Γ,uh[1])
# cell_u1 = get_cell_dof_values(uh[1])
# cell_id1 = _compute_cell_ids(uh[1],Γ)

# F2 = Gridap.MultiField.restrict_function(e,uh,2)
# g2 = _change_argument(hessian,F2,Γ,uh[2])
# cell_u2 = get_cell_dof_values(uh[2])
# cell_id2 = _compute_cell_ids(uh[2],Γ)

# agrad = i_to_y -> Gridap.Arrays.autodiff_array_gradient(g,i_to_y,cell_id)
# cell_grad = Gridap.Arrays.autodiff_array_jacobian(agrad,cell_u,cell_id)
# add_contribution!(terms,Γ,cell_grad)

# get_array(terms)[1]

# cell_u[1]
# cell_id

# get_cell_dof_values(uh)[2]

# j = hessian(x->e((x,uh[2])),uh[1])
# get_array(j)[1]

# j = hessian(x->e((dv[1],x)),uh[2])
# get_array(j)[1][1]

# j = hessian(e,uh)
# get_array(j)[1]

# function Gridap.Arrays.autodiff_array_hessian(a,i_to_x,j_to_i)
#   agrad = i_to_y -> Gridap.Arrays.autodiff_array_gradient(a,i_to_y,j_to_i)
#   Gridap.Arrays.autodiff_array_jacobian(agrad,i_to_x,j_to_i)
# end

# # dv = get_fe_basis(X)
# j_hes = hessian(e,uh)
# J_hes = assemble_matrix(j_hes,X,X)
