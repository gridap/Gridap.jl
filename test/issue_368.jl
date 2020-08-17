module TestIssue368
# Boundary integrals throw error trying to call zero(::Union{})
using Test
using LinearAlgebra
using Gridap

function solve_weak_dirichlet(model)
    V = TestFESpace(reffe=:Lagrangian, order=1, model=model, valuetype=Float64)

    u_e(pt) = -1/2*pt[1]^2

    F(u, v) = u ⊙ v - u_e ⊙ v
    B(u, v) = u ⊙ v - u_e ⊙ v

    trian = Triangulation(model)
    quad = CellQuadrature(trian, 2)
    t_inner = FETerm(F, trian, quad)

    btrian = BoundaryTriangulation(model, "boundary")
    bquad = CellQuadrature(btrian, 2)
    t_bdry = FETerm(B, btrian, bquad)

    U = TrialFESpace(V)
    op = FEOperator(U,V,t_inner, t_bdry)
    sol = solve(op)
    @test norm(integrate(sol - u_e, trian, quad)) < 1e3
    does_not_throw = true
end

model1d = CartesianDiscreteModel((0,1), 10)
model2d = CartesianDiscreteModel((0,1, 0, 1), (10,10))

@test_broken solve_weak_dirichlet(model1d)
@test solve_weak_dirichlet(model2d)
end
