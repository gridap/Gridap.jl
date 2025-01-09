
using Gridap
using Gridap.Geometry, Gridap.CellData

using Gridap.Geometry: HybridTriangulation, get_faces
using Gridap.CellData: change_domain

model = CartesianDiscreteModel((0,1,0,1),(4,4))

Ω = Triangulation(model)
Γ = Triangulation(ReferenceFE{1},model)

cell_to_bgcell = [6,7,10,11,14,15]
Λ = HybridTriangulation(model,cell_to_bgcell)
face_to_bgface = Λ.glue.face_glue.face_to_bgface


Ωv = view(Ω,[1,2,3,4])
is_change_possible(Ω,Ωv)
is_change_possible(Ωv,Ω)
num_cells(best_target(Ω,Ωv))

∂Ω = Boundary(model)
is_change_possible(Ω,∂Ω)
is_change_possible(∂Ω,Ω)


reffe = ReferenceFE(lagrangian,Float64,1)
VΩ = FESpace(Ω,reffe)
VΓ = FESpace(Γ,reffe)
X = MultiFieldFESpace([VΩ,VΓ])

u1 = zero(VΩ)
u2 = zero(VΓ)

u1*u2


n = get_normal_vector(Λ)

ΠΛ(u) = change_domain(u,Λ.face_trian,DomainStyle(u))
dΛ = Measure(Λ,2)
mass(u,v) = ∫(u*v)*dΛ
function mf_mass((_u1,_u2),(_v1,_v2))
  u1, u2, v1, v2 = ΠΛ(_u1), ΠΛ(_u2), ΠΛ(_v1), ΠΛ(_v2)
  ∫(u1*v1 + u1*v2 + u2*v1 + u2*v2)*dΛ
end

assemble_matrix(mass,VΩ,VΩ)
assemble_matrix(mass,VΓ,VΓ)
assemble_matrix(mf_mass,X,X)


function CellData.Measure(trian::Geometry.HybridTriangulation,args...;kwargs...)
  CellData.CompositeMeasure(trian.cell_trian,trian.face_trian,args...;kwargs...)
end

u = get_trial_fe_basis(X)
v = get_fe_basis(X)

c = mf_mass(u,v)

for trian in get_domains(c)
  arr = get_contribution(c,trian)
  display(arr)
end

arr = get_contribution(c,Λ.cell_trian)

#


