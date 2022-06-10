module SkeletonCellFieldPairTests

using Test
using FillArrays
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
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

u(x) = sin(norm(x))
U = TrialFESpace(V)
uh = FEFunction(U,rand(num_free_dofs(U)))

cellu = get_cell_dof_values(uh)
cellu_dual = lazy_map(Gridap.Arrays.DualizeMap(ForwardDiff.gradient),cellu)
uh_dual = CellField(U,cellu_dual)
@test eltype(cellu_dual[1]) <: ForwardDiff.Dual

scfp = SkeletonCellFieldPair(uh,uh,Λ) # non-dualized version
scfp_dualplus = SkeletonCellFieldPair(uh_dual,uh,Λ)

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

end # module
