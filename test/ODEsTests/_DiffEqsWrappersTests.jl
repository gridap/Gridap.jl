module DiffEqsWrappersTests

using Test

using Gridap
using Gridap.ODEs

# using DifferentialEquations
# using Sundials
using Gridap.Algebra: NewtonRaphsonSolver
using Base.Iterators

# FE problem (heat eq) using Gridap
function fe_problem(u, n)

  f(t) = x -> ∂t(u)(x, t) - Δ(u(t))(x)

  domain = (0, 1, 0, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition)

  order = 1

  reffe = ReferenceFE(lagrangian, Float64, order)
  V0 = FESpace(
    model,
    reffe,
    conformity=:H1,
    dirichlet_tags="boundary",
  )
  U = TransientTrialFESpace(V0, u)

  Ω = Triangulation(model)
  degree = 2 * order
  dΩ = Measure(Ω, degree)

  a(u, v) = ∫(∇(v) ⋅ ∇(u))dΩ
  b(v, t) = ∫(v * f(t))dΩ
  m(u, v) = ∫(v * u)dΩ

  res(t, u, v) = a(u, v) + m(∂t(u), v) - b(v, t)
  jac(t, u, du, v) = a(du, v)
  jac_t(t, u, dut, v) = m(dut, v)

  op = TransientFEOperator(res, jac, jac_t, U, V0)

  U0 = U(0.0)
  uh0 = interpolate_everywhere(u(0.0), U0)
  u0 = get_free_dof_values(uh0)

  return op, u0

end

# Solving the heat equation using Gridap.ODEs and DiffEqs
tspan = (0.0, 1.0)

u(x, t) = t
u(t) = x -> u(x, t)

# ISSUE 1: When I choose n > 2, even though the problem that we will solve is
# linear, the Sundials solvers seems to have convergence issues in the nonlinear
# solver (?). Ut returns errors
# [IDAS ERROR]  IDACalcIC Newton/Linesearch algorithm failed to converge.
# ISSUE 2: When I pass `jac_prototype` the code gets stuck

n = 3 # cells per dim (2D)
op, u0 = fe_problem(u, n)

# Some checks
res!, jac!, mass!, stif! = diffeq_wrappers(op)
J = prototype_jacobian(op, u0)
r = copy(u0)
θ = 1.0
t0 = 0.0
tF = 1.0
dt = 0.1
tθ = 1.0
dtθ = dt * θ
res!(r, u0, u0, nothing, tθ)
jac!(J, u0, u0, nothing, (1 / dtθ), tθ)

K = prototype_jacobian(op, u0)
M = prototype_jacobian(op, u0)
stif!(K, u0, u0, nothing, tθ)
mass!(M, u0, u0, nothing, tθ)
# Here you have the mass matrix M

@test (1 / dtθ) * M + K ≈ J

# To explore the Sundials solver options, e.g., BE with fixed time step dtd
f_iip = DAEFunction{true}(res!; jac=jac!)#, jac_prototype=J)
# jac_prototype is the way to pass my pre-allocated jacobian matrix
prob_iip = DAEProblem{true}(f_iip, u0, u0, tspan, differential_vars=[true, true, true, true])
# When I pass `jac_prototype` the code get stuck here:
# sol_iip = Sundials.solve(prob_iip, IDA(), reltol = 1e-8, abstol = 1e-8)
# @show sol_iip.u

# or iterator version
# integ = init(prob_iip, IDA(), reltol = 1e-8, abstol = 1e-8)
# step!(integ)

# Show using integrators as iterators
# for i in take(integ, 100)
# @show integ.u
# end

end # module DiffEqsWrappersTests

# FUTURE WORK: Check other options, not only Sundials

# ISSUE: Future work, add own (non)linear solver.
# Let us assume that we just want to consider a fixed point algorithm
# and we consider an implicit time integration of a nonlinear PDE.
# Our solvers are efficient since they re-use cache among
# nonlinear steps (e.g., symbolic factorization, arrays for numerical factorization)
# and for the transient case in our library, we can also reuse all this
# between time steps. Could we attain something like this using
# DifferentialEquations/Sundials?
# @ChrisRackauckas suggests to take a look at:
# https://docs.sciml.ai/latest/tutorials/advanced_ode_example/ shows swapping out linear solvers.
# https://docs.sciml.ai/latest/features/linear_nonlinear/ is all of the extra details.

# Try to pass solver too
# ls = LUSolver()
# ls_cache = nothing
# x = copy(u0)
# solve!(x,J,r,ls_cache)
#
