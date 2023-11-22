module TransientFESpacesTests

using Test

using Gridap
using Gridap.ODEs

import Gridap: ∇
import Gridap.ODEs: ∂t
import Gridap.ODEs: ∂tt

u(x, t) = (x[1] + x[2]) * t
u(t::Real) = x -> u(x, t)

∇u(x, t) = VectorValue(1, 1) * t
∇u(t::Real) = x -> ∇u(x, t)
∇(::typeof(u)) = ∇u

∂tu(t) = x -> x[1] + x[2]
∂t(::typeof(u)) = ∂tu

∂ttu(t) = x -> zero(x[1])
∂tt(::typeof(u)) = ∂ttu

# f(t) = x -> (x[1] + x[2])

domain = (0, 1, 0, 1)
partition = (5, 5)
model = CartesianDiscreteModel(domain, partition)

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
U = TransientTrialFESpace(V, u)
@test test_transient_trial_fe_space(U)

Ut = ∂t(U)
Utt = ∂tt(U)

for t in (0.0, 1.0, 2.0)
  # Dirichlet values of U
  _U0 = TrialFESpace(V, u(t))
  U0 = U(t)

  _ud0 = get_dirichlet_dof_values(_U0)
  ud0 = get_dirichlet_dof_values(U0)

  @test all(ud0 .≈ _ud0)

  # Dirichlet values of ∂t(U)
  _Ut0 = TrialFESpace(V, ∂tu(t))
  Ut0 = Ut(t)

  _utd0 = get_dirichlet_dof_values(_Ut0)
  utd0 = get_dirichlet_dof_values(Ut0)

  @test all(utd0 .≈ _utd0)

  # Dirichlet values of ∂tt(U)
  _Utt0 = TrialFESpace(V, ∂ttu(t))
  Utt0 = Utt(t)

  _uttd0 = get_dirichlet_dof_values(_Utt0)
  uttd0 = get_dirichlet_dof_values(Utt0)

  @test all(uttd0 .≈ _uttd0)
end

end # module TransientFESpacesTests
