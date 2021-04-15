module BoundaryDiscreteModelsTests

using Test
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Arrays
using Gridap.TensorValues

domain = (0,1,0,1,0,1)
cells = (3,3,3)
bgmodel = CartesianDiscreteModel(domain,cells)

labels = get_face_labeling(bgmodel)
bgface_to_mask = get_face_mask(labels,[22,23],2)

model = BoundaryDiscreteModel(Polytope{2},bgmodel,bgface_to_mask)

trian = Triangulation(model)
btrian = BoundaryTriangulation(model)
strian = SkeletonTriangulation(model)

@test get_background_triangulation(trian) === Triangulation(bgmodel)
@test get_background_triangulation(btrian) === trian
@test get_background_triangulation(strian) === trian
@test get_background_triangulation(Triangulation(model.model)) === Triangulation(model.model)

reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(model,reffe)
bgV = FESpace(bgmodel,reffe)

vh = FEFunction(V,rand(num_free_dofs(V)))
bgvh = FEFunction(bgV,rand(num_free_dofs(bgV)))

@test get_triangulation(vh) === Triangulation(model)
@test get_background_triangulation(get_triangulation(vh)) === Triangulation(bgmodel)
@test get_background_triangulation(get_triangulation(vh)) === get_triangulation(bgvh)

Ω = Triangulation(bgmodel)
Γ = trian
∂Γ = btrian
sΓ = strian

dΩ = Measure(Ω,3)
dΓ = Measure(Γ,3)
d∂Γ = Measure(∂Γ,3)
dsΓ = Measure(sΓ,3)

n_Γ = get_normal_vector(Γ)
n_∂Γ = get_normal_vector(∂Γ)

@test ∑(∫(1)dΓ) ≈ 2.0
@test ∑(∫(1)d∂Γ) ≈ 6.0
@test ∑(∫(1)dsΓ) ≈ 9.0

@test ∑(∫(n_Γ⋅VectorValue(1.0,0.0,0.0))dΓ) ≈ 0.0
@test ∑(∫(n_Γ⋅VectorValue(0.0,1.0,0.0))dΓ) ≈ -1.0
@test ∑(∫(n_Γ⋅VectorValue(0.0,0.0,1.0))dΓ) ≈ 1.0
@test ∑(∫(n_Γ⋅VectorValue(0.0,0.0,1.0))d∂Γ) ≈ 3.0
@test ∑(∫(n_∂Γ⋅VectorValue(0.0,0.0,1.0))d∂Γ) ≈ -1.0

end
