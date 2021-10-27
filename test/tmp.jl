using Pkg
Pkg.activate(".")

module tmp

using Gridap
using Gridap.Geometry

domain = (0, 1, 0, 1)
partition = (4,4)
model_Ω = CartesianDiscreteModel(domain,partition)# |> simplexify
labels_Ω = get_face_labeling(model_Ω)
add_tag_from_tags!(labels_Ω,"surface",[6])
Ω = Interior(model_Ω)
Γ = Boundary(model_Ω,tags="surface")
# writevtk(Ω,"O")
# writevtk(Γ,"G")

# function is_plate(x)
#   is_in = ([0.25 <= xm[1] <= 0.75  for xm in x])
#   minimum(is_in)
# end
# xΓ = get_cell_coordinates(Γ)
# Γb_to_Γ_mask = lazy_map(is_plate,xΓ)
# Γb_to_Γ = findall(Γb_to_Γ_mask)
# Γf_to_Γ = findall(!,Γb_to_Γ_mask)
# Γb = Triangulation(Γ,Γb_to_Γ)
# Γf = Triangulation(Γ,Γf_to_Γ)
# writevtk(Γb,"Gb")
# writevtk(Γf,"Gf")

# Λb = Skeleton(Γb)
# writevtk(Λb,"Lb")

# # Measures
# order = 2
# degree = 2*order
# dΩ = Measure(Ω,degree)
# dΓb = Measure(Γb,degree)
# dΓf = Measure(Γf,degree)
# dΛb = Measure(Λb,degree)

# # Normals
# nΛb = get_normal_vector(Λb)

# # FE spaces
# reffe = ReferenceFE(lagrangian,Float64,order)
# V_Ω = TestFESpace(Ω, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
# V_Γ = TestFESpace(Γ, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
# U_Ω = TrialFESpace(V_Ω)
# U_Γ = TrialFESpace(V_Γ)
# X = MultiFieldFESpace([U_Ω,U_Γ])
# Y = MultiFieldFESpace([V_Ω,V_Γ])

# # Weak form
# a((ϕ,η),(w,v)) = ∫( ∇(ϕ)⋅∇(w) )dΩ +
#   ∫( (η-im*ϕ)*(v + w) + im*η*w )dΓf + ∫( (η-im*ϕ)*v + im*η*w + Δ(η)*Δ(v) )dΓb +
#   ∫( jump(∇(v)⋅nΛb)*mean(Δ(η)) - mean(Δ(v))*jump(∇(η)⋅nΛb) + jump(∇(v)⋅nΛb)*jump(∇(η)⋅nΛb) )dΛb
# b((w,v)) =  0.0
# op = AffineFEOperator(a,b,X,Y)

end
