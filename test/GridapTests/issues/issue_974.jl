using Gridap

T = ComplexF64

# Operator creating a linear system where the matrix and RHS are Float64 but the unknown vector is ComplexF64.
domain = (0, 1, 0, 1)
partition = (4, 4)
model = CartesianDiscreteModel(domain, partition)

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)

V1 = TestFESpace(model, reffe; conformity=:H1, vector_type=Vector{T})
V2 = TestFESpace(model, reffe; conformity=:H1, vector_type=Vector{T})

U1 = TrialFESpace(V1)
U2 = TrialFESpace(V2)

Y = MultiFieldFESpace([V1, V2])
X = MultiFieldFESpace([U1, U2])

degree = 2 * order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

# Project constant 1 into both spaces.
a((u1, u2), (v1, v2)) = ∫(v1 * u1)dΩ + ∫(v2 * u2)dΩ
l((v1, v2)) = ∫(v1 + v2)dΩ

op = AffineFEOperator(a, l, X, Y)
uh1, uh2 = solve(op)
