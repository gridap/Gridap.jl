module HDGTests

using Test
using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.Fields, Gridap.Helpers

# u(x) = sin(2*π*x[1])*sin(2*π*x[2])*(1-x[1])*x[2]*(1-x[2])
u(x) = x[1] + x[2]
q(x) = -∇(u)(x)
f(x) = (∇ ⋅ q)(x)

D = 3
nc = Tuple(fill(4, D))
domain = Tuple(repeat([0, 1], D))
model = simplexify(CartesianDiscreteModel(domain,nc))

Ω = Triangulation(ReferenceFE{D}, model)
Γ = Triangulation(ReferenceFE{D-1}, model)

ptopo = Geometry.PatchTopology(model)
Ωp = Geometry.PatchTriangulation(model,ptopo)
Γp = Geometry.PatchBoundaryTriangulation(model,ptopo)

# Reference FEs
order = 1
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
X = MultiFieldFESpace([V, Q, M];style=mfs)

degree = 2*(order+1)
dΩp = Measure(Ωp,degree)
dΓp = Measure(Γp,degree)

τ = 1.0 # HDG stab parameter
n = get_normal_vector(Γp)
Πn(u) = u⋅n
Π(u) = change_domain(u,Γp,DomainStyle(u))
a((qh,uh,sh),(vh,wh,lh)) = ∫( qh⋅vh - uh*(∇⋅vh) - qh⋅∇(wh) )dΩp + ∫(sh*Πn(vh))dΓp +
                           ∫((Πn(qh) + τ*(Π(uh) - sh))*(Π(wh) + lh))dΓp
l((vh,wh,lh)) = ∫( f*wh )*dΩp

op = MultiField.StaticCondensationOperator(ptopo,X,a,l)
qh, uh, sh = solve(op)

dΩ = Measure(Ω,degree)
eh = uh - u
l2_uh = sqrt(sum(∫(eh⋅eh)*dΩ))
@test l2_uh < 1e-10

end # module