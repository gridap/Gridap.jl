module SkeletonCellFieldPairTests

using Test
using FillArrays
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.MultiField
using ForwardDiff

model = CartesianDiscreteModel((0.,1.,0.,1.),(3,3))
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)

dΩ = Measure(Ω,2)
dΓ = Measure(Γ,2)
dΛ = Measure(Λ,2)

n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)

order = 2
u(x) = sin(norm(x))
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe,conformity=:L2)
U = TrialFESpace(V)
uh = FEFunction(U,rand(num_free_dofs(U)))

cellu = get_cell_dof_values(uh)
cellu_dual = lazy_map(DualizeMap(ForwardDiff.gradient),cellu)
uh_dual = CellField(U,cellu_dual)
@test eltype(cellu_dual[1]) <: ForwardDiff.Dual

scfp = SkeletonCellFieldPair(uh,uh) # non-dualized version
scfp_dualplus = SkeletonCellFieldPair(uh_dual,uh)

# construction tests
@test get_triangulation(scfp) === get_triangulation(uh)
@test scfp.plus === uh.plus
@test scfp.minus === uh.minus

@test get_triangulation(scfp_dualplus) === get_triangulation(uh)
@test scfp_dualplus.plus === uh_dual.plus
@test scfp_dualplus.minus === uh.minus

# few basic tests if integrations on trians are done fine
@test sum(∫(scfp)*dΓ) == sum(∫(uh)*dΓ)
@test sum(∫(scfp)*dΩ) == sum(∫(uh)*dΩ)
@test sum(∫(mean(scfp))*dΛ) == sum(∫(mean(uh))*dΛ)
@test sum(∫(jump(scfp*n_Λ))*dΛ) == sum(∫(jump(uh*n_Λ))*dΛ)
@test sum(∫(uh*uh)*dΩ) == sum(∫(scfp*scfp)*dΩ)
@test sum(∫(u*uh)*dΩ) == sum(∫(u*scfp)*dΩ)
@test sum(∫(u-uh)*dΩ) == sum(∫(u-scfp)*dΩ)

# integrations involving SkeletonCellFieldPair derivative terms
@test sum(∫(∇(scfp))*dΩ) == sum(∫(∇(uh))*dΩ)
@test sum(∫(∇(scfp)⋅n_Γ)*dΓ) == sum(∫(∇(uh)⋅n_Γ)*dΓ)
@test sum(∫(mean(∇(uh)))*dΛ) == sum(∫(mean(∇(scfp)))*dΛ)
@test sum(∫(jump(∇(scfp)))*dΛ) == sum(∫(jump(∇(uh)))*dΛ)

# dualized version tests

# few basic tests if integrations on trians are done fine
@test sum(∫(scfp_dualplus)*dΓ).value == sum(∫(uh)*dΓ)
@test sum(∫(scfp_dualplus)*dΓ).partials == sum(∫(uh_dual)*dΓ).partials
@test sum(∫(scfp_dualplus)*dΩ).value == sum(∫(uh)*dΩ)
@test sum(∫(scfp_dualplus)*dΩ).partials == sum(∫(uh_dual)*dΩ).partials
@test sum(∫(mean(scfp_dualplus))*dΛ).value == sum(∫(mean(uh))*dΛ)
@test sum(∫(jump(scfp_dualplus*n_Λ))*dΛ)[1].value == sum(∫(jump(uh*n_Λ))*dΛ)[1]
@test sum(∫(uh*uh)*dΩ) == sum(∫(scfp_dualplus*scfp_dualplus)*dΩ).value
@test sum(∫(uh_dual*uh_dual)*dΩ).partials == sum(∫(scfp_dualplus*scfp_dualplus)*dΩ).partials
@test sum(∫(u*uh)*dΩ) == sum(∫(u*scfp_dualplus)*dΩ).value
@test sum(∫(u*uh_dual)*dΩ).partials == sum(∫(u*scfp_dualplus)*dΩ).partials
@test sum(∫(u-uh)*dΩ) == sum(∫(u-scfp_dualplus)*dΩ).value
@test sum(∫(u-uh_dual)*dΩ).partials == sum(∫(u-scfp_dualplus)*dΩ).partials

# integrations involving SkeletonCellFieldPair derivative terms
@test sum(∫(∇(scfp_dualplus))*dΩ)[1].value == sum(∫(∇(uh))*dΩ)[1]
@test sum(∫(∇(scfp_dualplus))*dΩ)[1].partials == sum(∫(∇(uh_dual))*dΩ)[1].partials
@test sum(∫(∇(scfp_dualplus)⋅n_Γ)*dΓ).value == sum(∫(∇(uh)⋅n_Γ)*dΓ)
@test sum(∫(∇(scfp_dualplus)⋅n_Γ)*dΓ).partials == sum(∫(∇(uh_dual)⋅n_Γ)*dΓ).partials
@test sum(∫(mean(∇(uh)))*dΛ)[1] == sum(∫(mean(∇(scfp_dualplus)))*dΛ)[1].value
@test sum(∫(jump(∇(scfp_dualplus)))*dΛ)[1].value == sum(∫(jump(∇(uh)))*dΛ)[1]

# AD of Integrals over Skeleton faces tested in FESpaces -> FEAutodiffTests.jl

# testing the MultiField SkeletonCellFieldPair (SCFP)
# this uses some low-level API to construct MultiFieldCellFields with SCFP

Q = TestFESpace(model,reffe,conformity=:L2)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

xh = FEFunction(X,rand(num_free_dofs(X)))
uh,ph = xh

cellu = lazy_map(DensifyInnerMostBlockLevelMap(),get_cell_dof_values(xh))
cellu_dual = lazy_map(DualizeMap(ForwardDiff.gradient),cellu)
single_fields_plus = SkeletonCellFieldPair[]
single_fields = SkeletonCellFieldPair[]
single_fields_dual = GenericCellField[]
U = xh.fe_space
nfields = length(U.spaces)
cell_dofs_field_offsets = MultiField._get_cell_dofs_field_offsets(xh)

for i in 1:nfields
  view_range = cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
  cell_values_field = lazy_map(a->view(a,view_range),cellu_dual)
  cf_dual = CellField(U.spaces[i],cell_values_field)
  scfp_plus = SkeletonCellFieldPair(cf_dual, xh[i])
  scfp = SkeletonCellFieldPair(xh[i], xh[i])
  push!(single_fields_dual,cf_dual)
  push!(single_fields_plus,scfp_plus)
  push!(single_fields,scfp)
end

xh_scfp = MultiFieldCellField(single_fields)
uh_scfp, ph_scfp = xh_scfp

xh_scfp_dual_plus = MultiFieldCellField(single_fields_plus)
uh_scfp_dual_plus, ph_scfp_dual_plus = xh_scfp_dual_plus

xh_dual = MultiFieldCellField(single_fields_dual)
uh_dual, ph_dual = xh_dual

g_Λ((uh,ph)) = ∫( mean(uh) + mean(ph) + mean(uh)*mean(ph) )dΛ
a_Λ((uh,ph)) = ∫( - jump(uh*n_Λ)⊙mean(∇(ph))
                  - mean(∇(uh))⊙jump(ph*n_Λ)
                  + jump(uh*n_Λ)⊙jump(ph*n_Λ) )dΛ

# few basic tests if integrations on trians are done fine
@test sum(∫(uh_scfp*ph_scfp)*dΓ) == sum(∫(uh*ph)*dΓ)
@test sum(∫(uh_scfp - ph_scfp)*dΓ) == sum(∫(uh - ph)*dΓ)
@test sum(∫(uh_scfp*ph_scfp)*dΩ) == sum(∫(uh*ph)*dΩ)
@test sum(∫(uh_scfp - ph_scfp)*dΩ) == sum(∫(uh - ph)*dΩ)

# integrations involving SkeletonCellFieldPair derivative terms
@test sum(∫(ph_scfp*∇(uh_scfp))*dΩ) == sum(∫(ph*∇(uh))*dΩ)
@test sum(∫(uh_scfp*∇(ph_scfp))*dΩ) == sum(∫(uh*∇(ph))*dΩ)
@test sum(∫(uh_scfp*∇(ph_scfp)⋅n_Γ)*dΓ) == sum(∫(uh*∇(ph)⋅n_Γ)*dΓ)
@test sum(∫(ph_scfp*∇(uh_scfp)⋅n_Γ)*dΓ) == sum(∫(ph*∇(uh)⋅n_Γ)*dΓ)

# dualized tests

# few basic tests if integrations on trians are done fine
@test sum(∫(uh_scfp_dual_plus*ph_scfp_dual_plus)*dΓ).value == sum(∫(uh*ph)*dΓ)
@test sum(∫(uh_scfp_dual_plus*ph_scfp_dual_plus)*dΓ).partials == sum(∫(uh_dual*ph_dual)*dΓ).partials
@test sum(∫(uh_scfp_dual_plus - ph_scfp_dual_plus)*dΓ).value == sum(∫(uh - ph)*dΓ)
@test sum(∫(uh_scfp_dual_plus - ph_scfp_dual_plus)*dΓ).partials == sum(∫(uh_dual - ph_dual)*dΓ).partials
@test sum(∫(uh_scfp_dual_plus*ph_scfp_dual_plus)*dΩ).value == sum(∫(uh*ph)*dΩ)
@test sum(∫(uh_scfp_dual_plus*ph_scfp_dual_plus)*dΩ).partials == sum(∫(uh_dual*ph_dual)*dΩ).partials
@test sum(∫(uh_scfp_dual_plus - ph_scfp_dual_plus)*dΩ).value == sum(∫(uh - ph)*dΩ)
@test sum(∫(uh_scfp_dual_plus - ph_scfp_dual_plus)*dΩ).partials == sum(∫(uh_dual - ph_dual)*dΩ).partials

# integrations involving SkeletonCellFieldPair derivative terms
@test sum(∫(ph_scfp_dual_plus*∇(uh_scfp_dual_plus))*dΩ)[1].value == sum(∫(ph*∇(uh))*dΩ)[1]
@test sum(∫(ph_scfp_dual_plus*∇(uh_scfp_dual_plus))*dΩ)[1].partials == sum(∫(ph_dual*∇(uh_dual))*dΩ)[1].partials
@test sum(∫(uh_scfp_dual_plus*∇(ph_scfp_dual_plus))*dΩ)[1].value == sum(∫(uh*∇(ph))*dΩ)[1]
@test sum(∫(uh_scfp_dual_plus*∇(ph_scfp_dual_plus))*dΩ)[1].partials == sum(∫(uh_dual*∇(ph_dual))*dΩ)[1].partials
@test sum(∫(uh_scfp_dual_plus*∇(ph_scfp_dual_plus)⋅n_Γ)*dΓ).value == sum(∫(uh*∇(ph)⋅n_Γ)*dΓ)
@test sum(∫(uh_scfp_dual_plus*∇(ph_scfp_dual_plus)⋅n_Γ)*dΓ).partials == sum(∫(uh_dual*∇(ph_dual)⋅n_Γ)*dΓ).partials
@test sum(∫(ph_scfp_dual_plus*∇(uh_scfp_dual_plus)⋅n_Γ)*dΓ).value == sum(∫(ph*∇(uh)⋅n_Γ)*dΓ)
@test sum(∫(ph_scfp_dual_plus*∇(uh_scfp_dual_plus)⋅n_Γ)*dΓ).partials == sum(∫(ph_dual*∇(uh_dual)⋅n_Γ)*dΓ).partials

# test for SkeletonTriangulation
@test sum(g_Λ((uh_scfp_dual_plus,ph_scfp_dual_plus))).value == sum(g_Λ((uh,ph)))
@test sum(a_Λ((uh_scfp_dual_plus,ph_scfp_dual_plus))).value == sum(a_Λ((uh,ph)))

end # module
