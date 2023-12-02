module ODESolversTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")
include("ODESolversMocks.jl")

t0 = randn()
tF = t0 + rand()
dt = (tF - t0) / 10

num_eqs = 2

M = randn(num_eqs, num_eqs)
C = randn(num_eqs, num_eqs)
K = randn(num_eqs, num_eqs)
while iszero(det(M + dt * K)) || iszero(det(M + dt * C + 3 * dt^2 / 2 * K))
  M = randn(num_eqs, num_eqs)
  C = randn(num_eqs, num_eqs)
  K = randn(num_eqs, num_eqs)
end
α = randn(num_eqs)
f(t) = -exp.(α .* t)

u0 = randn(num_eqs)
v0 = randn(num_eqs)

_odeop1 = ODEOperatorMock1
_odeslvr1 = ODESolverMock1
_disop1 = DiscreteODEOperatorMock1
forms1 = (M, K)
us01 = (u0,)
_r1(x) = M * x + K * (u0 + dt * x) + f(t0 + dt)
_J1(x) = M + dt * K
x̃1 = -(M + dt * K) \ (K * u0 + f(t0 + dt))

_odeop2 = ODEOperatorMock2
_odeslvr2 = ODESolverMock2
_disop2 = DiscreteODEOperatorMock2
forms2 = (M, C, K)
us02 = (u0, v0,)
_r2(x) = M * x + C * (v0 + dt * x) + K * (u0 + dt * v0 + 3 * dt^2 / 2 * x) + f(t0 + dt)
_J2(x) = M + dt * C + 3 * dt^2 / 2 * K
x̃2 = -(M + dt * C + 3 * dt^2 / 2 * K) \ (C * v0 + K * (u0 + dt * v0) + f(t0 + dt))

disslvr = NLSolverMock()

for (_odeop, _odeslvr, _disop, forms, us0, _r, _J, x̃) in (
  (_odeop1, _odeslvr1, _disop1, forms1, us01, _r1, _J1, x̃1),
  (_odeop2, _odeslvr2, _disop2, forms2, us02, _r2, _J2, x̃2)
)
  for T in (NonlinearODE, MassLinearODE, LinearODE)
    odeop = _odeop{T}(forms..., f)
    disop = _disop(odeop, nothing, t0, us0, dt)

    # MockDiscreteODEOperator tests
    x = randn(num_eqs)
    r̃ = _r(x)
    J̃ = _J(x)

    r = allocate_residual(disop, x)
    residual!(r, disop, x)
    @test r ≈ r̃

    J = allocate_jacobian(disop, x)
    fill!(J, 0)
    jacobian!(J, disop, x)
    @test J ≈ J̃

    # NLSolverMock tests
    disslvrcache = solve!(x, disslvr, disop)
    r, J, du = disslvrcache
    @test r ≈ r̃
    @test J ≈ J̃
    @test x ≈ x̃

    # ODESolver tests
    odeslvr = _odeslvr(disslvr, dt)

    # The following is a little ugly because the input of a first-order ODE
    # solver is a vector and not a tuple with one element
    usF = copy.(us0)
    usF, tF, cache = solve_step!(usF, odeslvr, odeop, us0, t0, nothing)

    if length(us0) == 1
      u0, = us0
      uF, = usF
      @test tF ≈ t0 + dt
      @test uF ≈ u0 + dt * x̃
    elseif length(us0) == 2
      u0, v0 = us0
      uF, vF = usF
      @test tF ≈ t0 + dt
      @test vF ≈ v0 + dt * x̃
      @test uF ≈ u0 + dt * v0 + 3 * dt^2 / 2 * x̃
    end

    @test test_ode_solver(odeslvr, odeop, t0, us0, dt)
  end
end

end # module ODESolversTests
