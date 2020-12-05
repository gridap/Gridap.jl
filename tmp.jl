using Gridap
using Gridap.FESpaces

domain = (0,2,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"inlet",[7])
add_tag_from_tags!(labels,"outlet",[8])
add_tag_from_tags!(labels,"free_surface",[3,4,6])
add_tag_from_tags!(labels,"water",[9])
writevtk(model,"tmp")

model_1d = DiscreteModel(Polytope{1},model)
model_fs = DiscreteModel(model_1d,tags="free_surface")

reffe =  ReferenceFE(:Lagrangian,Float64,1)
V = TestFESpace(
        model,
        reffe,
        conformity = :H1
)
U = TrialFESpace(V)
V_η = TestFESpace(
        model,
        reffe,
        conformity = :H1,
        dirichlet_tags=["water","inlet","bottom","outlet"]
)
U_η = TrialFESpace(V_η,0.0)

Y = MultiFieldFESpace([V,V_η])
X = MultiFieldFESpace([U,U_η])

degree = 2
Ω = Triangulation(model)
dΩ = LebesgueMeasure(Ω, degree)

# Neumann BCs
Ninlet = "inlet"
Nbottom = "bottom"
Noutlet = "outlet"
FS = "free_surface"

Γᵢ = BoundaryTriangulation(model, tags=Ninlet)
dΓᵢ = LebesgueMeasure(Γᵢ, degree)

Γb = BoundaryTriangulation(model, tags=Nbottom)
dΓb = LebesgueMeasure(Γb, degree)

Γₒ = BoundaryTriangulation(model, tags=Noutlet)
dΓₒ = LebesgueMeasure(Γₒ, degree)

Γfs = BoundaryTriangulation(model, tags=FS)
dΓfs = LebesgueMeasure(Γfs, degree)

U0 = 2.0
hInlet(x) = U0
hBottom(x) = 0
hOutlet(x) = 0

a_Ω(ϕ,v) = ∫( ∇(ϕ) ⊙ ∇(v) )dΩ
a_Γ((ϕ,η),(v,w)) = ∫( -(w*η) + 0.5*(v*ϕ + v*9.81*η + w/9.81*ϕ + w*η) )dΓfs
l_Γ(v) = ∫( v * hInlet )dΓᵢ + ∫( v * hBottom )dΓb + ∫( v * hOutlet )dΓₒ

a((ϕ,η),(v,w)) = a_Ω(ϕ,v) + a_Γ((ϕ,η),(v,w))
l((v,w)) = l_Γ(v)

opϕη = AffineFEOperator(a, l, X, Y)
sol = solve(opϕη)

writevtk(Ω,"restmp",cellfields=["phi"=>sol[1],"eta"=>sol[2]])
