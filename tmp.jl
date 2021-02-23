using Gridap
using Gridap.FESpaces
using Gridap.Geometry

domain = (0,2,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

bgface_to_mask = get_face_mask(labels,[3,4,6],1)
model_fs =BoundaryDiscreteModel(Polytope{1},model,bgface_to_mask)

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

Γfs = Triangulation(model_fs)
dΓfs = Measure(Γfs, degree)

U0 = 2.0
hInlet(x) = U0
hBottom(x) = 0
hOutlet(x) = 0

a_Ω(ϕ,w) = ∫( ∇(ϕ) ⊙ ∇(w) )dΩ
#a_Γ((ϕ,η),(v,w)) = ∫( -(w*η) + 0.5*(v*ϕ + v*9.81*η + w/9.81*ϕ + w*η) )dΓfs
a_Γ((ϕ,η),(w,v)) = ∫( v*η )dΓfs #+ w*η )dΓfs

a((ϕ,η),(w,v)) = a_Ω(ϕ,w) + a_Γ((ϕ,η),(w,v))
l((w,v)) = 0.0

opϕη = AffineFEOperator(a, l, X, Y)
sol = solve(opϕη)

writevtk(Ω,"restmp",cellfields=["phi"=>sol[1],"eta"=>sol[2]])
