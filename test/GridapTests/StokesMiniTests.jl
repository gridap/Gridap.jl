module StokesMiniTests

using Test
using Gridap

# Make sure `∇(b)` and `Δ(b)` are both continuous.
function b((x, y))
	x0 = floor(x)
	y0 = floor(y)
	ret = (x - x0) * (y - y0) * (x0 + 1 - x) * (y0 + 1 - y)
	sign = (-1)^(x0+y0)
	return 6^4 * sign * ret
end

u(x) = VectorValue(b(x) + 2*x[2], 3*b(x) + x[2] + 2*x[1])
p(x) = x[1] - 3*x[2]
f(x) = -Δ(u)(x) + ∇(p)(x)
g(x) = (∇⋅u)(x)
∇u(x) = ∇(u)(x)

domain = (0, 3, 0, 3)
partition = (3, 3)
model = CartesianDiscreteModel(domain, partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels, "dirichlet", [1, 2, 5])
add_tag_from_tags!(labels, "neumann", [3, 4, 6, 7, 8])

reffe_u = ReferenceFE(lagrangian, VectorValue{2, Float64}, 1)
reffe_b = ReferenceFE(bubble, VectorValue{2, Float64})
reffe_p = ReferenceFE(lagrangian, Float64, 1)

V = TestFESpace(model, reffe_u, labels = labels, dirichlet_tags = "dirichlet", conformity = :H1)
R = TestFESpace(model, reffe_b)
Q = TestFESpace(model, reffe_p, labels = labels, conformity = :H1, dirichlet_tags = "tag_1")

U = TrialFESpace(V, u)
B = TrialFESpace(R)
P = TrialFESpace(Q, p)

Y = MultiFieldFESpace([V, R, Q])
X = MultiFieldFESpace([U, B, P])

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model, labels, tags = "neumann")
n_Γ = get_normal_vector(Γ)

degree = 4
dΩ = Measure(Ω, degree)
dΓ = Measure(Γ, degree)

a((u, p), (v, q)) = ∫(∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u))*dΩ
l((v, q)) = ∫(v⋅f + q*g)*dΩ + ∫(v⋅(n_Γ⋅∇u) - (n_Γ⋅v)*p)*dΓ

mini_a((u, b, p), (v, r, q)) = a((u + b, p), (v+r, q))
mini_l((v, r, q)) = l((v+r, q))

op = AffineFEOperator(mini_a, mini_l, X, Y)

uch, ubh, ph = solve(op)
uh = ubh + uch

eu = u - uh
ep = p - ph

l2(u) = sqrt(sum(∫(u⊙u)*dΩ))
h1(u) = sqrt(sum(∫(u⊙u + ∇(u)⊙∇(u))*dΩ))

eu_l2 = l2(eu)
eu_h1 = h1(eu)
ep_l2 = l2(ep)

tol = 1.0e-10
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module