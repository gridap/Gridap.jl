using Gridap
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.Arrays

function is_beam(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  (0.5 <= x[1] <= 1.5 ) * ( x[2] ≈ 1 )
end

domain = (0,2,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)
Ω = Triangulation(model)

labels = get_face_labeling(model)
bgface_to_mask = get_face_mask(labels,[3,4,6],1)
Γface_to_bgface = findall(bgface_to_mask)
model_Γ = BoundaryDiscreteModel(Polytope{1},model,Γface_to_bgface)

Γ = Triangulation(model_Γ)
Γface_coords = get_cell_coordinates(Γ)
Γface_mask = lazy_map(is_beam,Γface_coords)
Γbface_Γface = findall(Γface_mask)
Γfface_Γface = findall(!,Γface_mask)
Γb = BoundaryTriangulation(model,view(Γface_to_bgface,Γbface_Γface))
Γf = BoundaryTriangulation(model,view(Γface_to_bgface,Γfface_Γface))

#writevtk(Ω,"Ω")
#writevtk(Γ,"Γ")
#writevtk(Γb,"Γb")
#writevtk(Γf,"Γf")

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V_ϕ = FESpace(model,reffe)
V_η = FESpace(model_Γ,reffe)
V = MultiFieldFESpace([V_ϕ,V_η])

degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
dΓb = Measure(Γb,degree)
dΓf = Measure(Γf,degree)

a_Ω(ϕ,w) = ∫( ∇(ϕ) ⊙ ∇(w) )dΩ
a_Γ((ϕ,η),(w,v)) =
  ∫( v*η + w*η + w*ϕ )dΓb +
  ∫( v*η + w*η + w*ϕ )dΓf +
  0*∫( v*η + w*η + w*ϕ )dΓ
a((ϕ,η),(w,v)) = a_Ω(ϕ,w) + a_Γ((ϕ,η),(w,v))
l((w,v)) = 0.0

op = AffineFEOperator(a,l,V,V)
ϕh,ηh = solve(op)

writevtk(Ω,"Ω",cellfields=["ϕh"=>ϕh])
writevtk(Γ,"Γ",cellfields=["ϕh"=>ϕh,"ηh"=>ηh])
writevtk(Γb,"Γb",cellfields=["ϕh"=>ϕh,"ηh"=>ηh])
writevtk(Γf,"Γf",cellfields=["ϕh"=>ϕh,"ηh"=>ηh])

