module ODEOperatorsTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")

t0 = randn()

num_eqs = 2

M = randn(num_eqs, num_eqs)
C = randn(num_eqs, num_eqs)
K = randn(num_eqs, num_eqs)
α = randn(num_eqs)
f(t) = -exp.(α .* t)

u0 = randn(num_eqs)
v0 = randn(num_eqs)
a0 = randn(num_eqs)

_odeop1 = ODEOperatorMock1
us1 = (u0, v0)
forms1 = (M, K)

_odeop2 = ODEOperatorMock2
us2 = (u0, v0, a0)
forms2 = (M, C, K)

for (_odeop, forms, us) in (
  (_odeop1, forms1, us1),
  (_odeop2, forms2, us2)
)
  for T in (NonlinearODE, MassLinearODE, LinearODE)
    odeop = _odeop{T}(forms..., f)

    odeopcache = allocate_odeopcache(odeop, t, us)
    update_odeopcache!(odeopcache, odeop, t)

    r = allocate_residual(odeop, t, us, odeopcache)
    @test size(r) == (num_eqs,)

    J = allocate_jacobian(odeop, t, us, odeopcache)
    @test size(J) == (num_eqs, num_eqs)

    _r = f(t)
    for (formi, ui) in zip(reverse(forms), us)
      _r .+= formi * ui
    end
    residual!(r, odeop, t, us, odeopcache)
    @test r ≈ _r

    _J = zeros(num_eqs, num_eqs)
    fill!(J, 0)
    for (i, formi) in enumerate(reverse(forms))
      _J .+= formi
      jacobian!(J, odeop, t, us, i - 1, 1, odeopcache)
      @test J ≈ _J
    end

    @test test_ode_operator(odeop, t, us)
  end
end

end # module ODEOperatorsTests
