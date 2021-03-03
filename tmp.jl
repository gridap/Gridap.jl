using Gridap
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.Arrays

domain = (0,2,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)
labels = get_face_labeling(model)

bgface_to_mask = get_face_mask(labels,[3,4,6],1)
model_fs =BoundaryDiscreteModel(Polytope{1},model,bgface_to_mask)

# Beam domain
function is_beam(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  (0.5 <= x[1] <= 1.5 ) * ( x[2] ≈ 1 )
end
∂Ω = BoundaryTriangulation(model,bgface_to_mask)
Γ_coords = get_cell_coordinates(∂Ω)
Γ_cell_to_is_beam = collect1d(lazy_map(is_beam,Γ_coords))
Γb_cell_to_Γ_cell = findall(Γ_cell_to_is_beam)
Γfs_cell_to_Γ_cell = findall(collect(Bool, .! Γ_cell_to_is_beam))
Γb = Triangulation(∂Ω,Γb_cell_to_Γ_cell)
Γfs = Triangulation(∂Ω,Γfs_cell_to_Γ_cell)

reffe =  ReferenceFE(lagrangian,Float64,1)
V = TestFESpace(
        model,
        reffe,
        conformity = :H1
)
U = TrialFESpace(V)
V_η = TestFESpace(
        model_fs,
        reffe,
        conformity = :H1
)
U_η = TrialFESpace(V_η)

Y = MultiFieldFESpace([V,V_η])
X = MultiFieldFESpace([U,U_η])

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

Γ = Triangulation(model_fs)
dΓ = Measure(Γ, degree)
dΓfs = Measure(Γfs, degree)
dΓb = Measure(Γb, degree)

U0 = 2.0
hInlet(x) = U0
hBottom(x) = 0
hOutlet(x) = 0

a_Ω(ϕ,w) = ∫( ∇(ϕ) ⊙ ∇(w) )dΩ
#a_Γ((ϕ,η),(v,w)) = ∫( -(w*η) + 0.5*(v*ϕ + v*9.81*η + w/9.81*ϕ + w*η) )dΓfs
#a_Γ((ϕ,η),(w,v)) = ∫( v*η )dΓfs #+ w*η )dΓfs
a_Γ((ϕ,η),(w,v)) = ∫( v*η + w*η + w*ϕ )dΓfs + ∫( v*η + w*η + w*ϕ )dΓb

a((ϕ,η),(w,v)) = a_Ω(ϕ,w) + a_Γ((ϕ,η),(w,v))
l((w,v)) = 0.0

opϕη = AffineFEOperator(a, l, X, Y)
sol = solve(opϕη)

writevtk(Ω,"restmp_vol",cellfields=["phi"=>sol[1]])
writevtk(Γfs,"restmp_surf",cellfields=["phi"=>sol[1],"eta"=>sol[2]])
