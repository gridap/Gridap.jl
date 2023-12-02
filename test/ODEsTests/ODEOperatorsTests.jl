module ODEOperatorsTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")

M = randn(2, 2)
C = randn(2, 2)
K = randn(2, 2)
f(t) = [cospi(t), sinpi(t)]

t = randn()
u = randn(2)
v = randn(2)
a = randn(2)

_odeop1 = ODEOperatorMock1
us1 = (u, v)
forms1 = (M, K)

_odeop2 = ODEOperatorMock2
us2 = (u, v, a)
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
    @test size(r) == (2,)

    J = allocate_jacobian(odeop, t, us, odeopcache)
    @test size(J) == (2, 2)

    _r = f(t)
    for (formi, ui) in zip(reverse(forms), us)
      _r .+= formi * ui
    end
    residual!(r, odeop, t, us, odeopcache)
    @test r ≈ _r

    _J = zeros(2, 2)
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
