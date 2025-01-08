
using Gridap
using Gridap.Geometry

using Gridap.Geometry: HybridTriangulation, get_faces
using Gridap.CellData: change_domain

model = CartesianDiscreteModel((0,1,0,1),(4,4))

Ω = Triangulation(model)
Γ = Triangulation(ReferenceFE{1},model)
Λ = HybridTriangulation(model,[6,7,10,11,14,15])

reffe = ReferenceFE(lagrangian,Float64,1)
VΩ = FESpace(Ω,reffe)
VΓ = FESpace(Γ,reffe)
X = MultiFieldFESpace([VΩ,VΓ])

n = get_normal_vector(Λ)

ΠΛ(u) = change_domain(u,Λ,DomainStyle(u))
dΛ = Measure(Λ,2)
mass(u,v) = ∫(u*v)*dΛ
function mf_mass((_u1,_u2),(_v1,_v2))
  u1, u2, v1, v2 = ΠΛ(_u1), ΠΛ(_u2), ΠΛ(_v1), ΠΛ(_v2)
  ∫(u1*v1 + u1*v2 + u2*v1 + u2*v2)*dΛ
end

assemble_matrix(mass,VΩ,VΩ)
assemble_matrix(mass,VΓ,VΓ)
assemble_matrix(mf_mass,X,X)
