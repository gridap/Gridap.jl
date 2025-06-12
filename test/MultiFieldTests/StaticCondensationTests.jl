module StaticCondensationTests

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers

u(x) = x[1]*(1-x[1])*x[2]*(1-x[2])
q(x) = -∇(u)(x)
f(x) = (∇ ⋅ q)(x)

nc = (8,8)
model = UnstructuredDiscreteModel(simplexify(CartesianDiscreteModel((0,1,0,1),nc)))
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

# Reference FEs
order = 2
reffeV = ReferenceFE(lagrangian, VectorValue{D, Float64}, order; space=:P)
reffeQ = ReferenceFE(lagrangian, Float64, order; space=:P)
reffeM = ReferenceFE(lagrangian, Float64, order; space=:P)

# HDG test FE Spaces
V_test = TestFESpace(Ω, reffeV; conformity=:L2)
Q_test = TestFESpace(Ω, reffeQ; conformity=:L2)
M_test = TestFESpace(Γ, reffeM; conformity=:L2, dirichlet_tags="boundary")

# HDG trial FE Spaces
V = TrialFESpace(V_test)
Q = TrialFESpace(Q_test)
M = TrialFESpace(M_test, u)

mfs = MultiField.BlockMultiFieldStyle(2,(2,1))
X_full = MultiFieldFESpace([V, Q, M];style=mfs)

τ = 1.0 # HDG stab parameter

degree = 2*(order+1)
dΩp = Measure(Ωp,degree)
dΓp = Measure(Γp,degree)

n = get_normal_vector(Γp)
Πn(u) = u⋅n
Π(u) = change_domain(u,Γp,DomainStyle(u))
a((qh,uh,sh),(vh,wh,lh)) = ∫( qh⋅vh - uh*(∇⋅vh) - qh⋅∇(wh) )dΩp + ∫(sh*Πn(vh))dΓp +
                           ∫((Πn(qh) + τ*(Π(uh) - sh))*(Π(wh) + lh))dΓp
l((vh,wh,lh)) = ∫( f*wh )*dΩp

op = MultiField.StaticCondensationOperator(ptopo,X_full,a,l)
qh, uh, sh = solve(op)

eu = uh - u
eq = qh - q
@test sum(∫(eu⋅eu)*dΩp) < 1.e-8
@test sum(∫(eq⋅eq)*dΩp) < 1.e-8

end # module