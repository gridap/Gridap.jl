module TestIssue349
using Test
using Gridap

# solve
# E = - ∇ϕ
# ρ = divergence(E)

ϕ_e(pt) = -1/2*pt[1]^2
E_e(pt) = pt
ρ(pt) = 1

model = CartesianDiscreteModel((-1,1), 10)
order = 1

V_E = TestFESpace(
  reffe=:Lagrangian, conformity=:H1, valuetype=VectorValue{1,Float64},
  model=model, order=order,
)

V_ϕ = TestFESpace(
  reffe=:Lagrangian, conformity=:H1, valuetype=Float64,
  model=model, order=order,  dirichlet_tags="boundary",
)

U_ϕ = TrialFESpace(V_ϕ, ϕ_e)
U_E = TrialFESpace(V_E)

V = MultiFieldFESpace([V_ϕ, V_E])
U = MultiFieldFESpace([U_ϕ, U_E])

function A_works(u, v)
    ϕ, E = u
    v_ϕ, v_E = v
    -∇(ϕ) ⊙ v_E - E ⊙ v_E + divergence(E) ⊙ v_ϕ
end

function A_fails(u, v)
    ϕ, E = u
    v_ϕ, v_E = v
    (-∇(ϕ) - E) ⊙ v_E + divergence(E) ⊙ v_ϕ
end

function b(v)
    v_ϕ, v_E = v
    ρ ⊙ v_ϕ
end

function errornorm(u_e, u, train, quad)
    arr = integrate((u_e - u)⊙(u_e - u), trian, quad)
    sqrt(sum(arr))
end

trian = Triangulation(model)
quad = CellQuadrature(trian, 2)
op_works = AffineFEOperator(U, V, AffineFETerm(A_works, b, trian, quad))
op_fails = AffineFEOperator(U, V, AffineFETerm(A_fails, b, trian, quad))
sol = solve(op_works) # works in #349
ϕ_sol, E_sol = sol
@test errornorm(ϕ_e, ϕ_sol, trian, quad) < 1e-2
@test errornorm(E_e, E_sol, trian, quad) < 1e-4
@test_broken begin
    ϕ_sol2, E_sol2 = solve(op_fails) # throws LinearAlgebra.SingularException(0) in #349
    @test errornorm(ϕ_e, ϕ_sol, trian, quad) ≈
        errornorm(ϕ_e, ϕ_sol2, trian, quad)
    @test errornorm(E_e, E_sol, trian, quad) ≈
        errornorm(E_e, E_sol2, trian, quad)
    true
end

end#module
