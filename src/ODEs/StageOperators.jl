##########################
# NonlinearStageOperator #
##########################
"""
    abstract type StageOperator <: NonlinearOperator end

Operator used to perform one stage within one time step of an `ODESolver`.

# Mandatory
- [`allocate_residual(nlop, x)`](@ref)
- [`residual!(r, nlop, x)`](@ref)
- [`allocate_jacobian(nlop, x)`](@ref)
- [`jacobian!(J, nlop, x)`](@ref)
"""
abstract type StageOperator <: NonlinearOperator end

##########################
# NonlinearStageOperator #
##########################
"""
    struct NonlinearStageOperator <: StageOperator end

Nonlinear stage operator representing `res(x) = residual(t, us(x)...) = 0`,
where `x` is the stage unknown and `us(x)` denotes the point where the residual
of the ODE is to be evaluated. It is assumed that the coordinates of `us(x)`
are linear in `x`, and the coefficients in front of `x` called `ws` are scalar,
i.e. `ws[k] = d/dx us[k](x)` is a scalar constant.
"""
struct NonlinearStageOperator <: StageOperator
  odeop::ODEOperator
  odeopcache
  tx::Real
  usx::Function
  ws::Tuple{Vararg{Real}}
end

# NonlinearOperator interface
function Algebra.allocate_residual(
  nlop::NonlinearStageOperator, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache
  tx = nlop.tx
  usx = nlop.usx(x)
  allocate_residual(odeop, tx, usx, odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  nlop::NonlinearStageOperator, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache
  tx = nlop.tx
  usx = nlop.usx(x)
  residual!(r, odeop, tx, usx, odeopcache)
end

function Algebra.allocate_jacobian(
  nlop::NonlinearStageOperator, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache
  tx = nlop.tx
  usx = nlop.usx(x)
  allocate_jacobian(odeop, tx, usx, odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  nlop::NonlinearStageOperator, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache
  tx = nlop.tx
  usx = nlop.usx(x)
  ws = nlop.ws
  jacobian!(J, odeop, tx, usx, ws, odeopcache)
  J
end

#######################
# LinearStageOperator #
#######################
"""
    struct LinearStageOperator <: StageOperator end

Linear stage operator representing `res(x) = J(t, us) x + r(t, us) = 0`,
where `x` is the stage unknown and `us` denotes the point where the residual
of the ODE is to be evaluated.
"""
struct LinearStageOperator <: StageOperator
  J::AbstractMatrix
  r::AbstractVector
  reuse::Bool
end

function LinearStageOperator(
  odeop::ODEOperator, odeopcache,
  tx::Real, usx::Tuple{Vararg{AbstractVector}},
  ws::Tuple{Vararg{Real}},
  J::AbstractMatrix, r::AbstractVector, reuse::Bool, sysslvrcache
)
  residual!(r, odeop, tx, usx, odeopcache)

  if isnothing(sysslvrcache) || !reuse
    jacobian!(J, odeop, tx, usx, ws, odeopcache)
  end

  LinearStageOperator(J, r, reuse)
end

# NonlinearOperator interface
function Algebra.allocate_residual(
  lop::LinearStageOperator, x::AbstractVector
)
  r = allocate_in_range(typeof(lop.r), lop.J)
  fill!(r, zero(eltype(r)))
  r
end

function Algebra.residual!(
  r::AbstractVector,
  lop::LinearStageOperator, x::AbstractVector
)
  mul!(r, lop.J, x)
  axpy!(1, lop.r, r)
  r
end

function Algebra.allocate_jacobian(
  lop::LinearStageOperator, x::AbstractVector
)
  lop.J
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  lop::LinearStageOperator, x::AbstractVector
)
  copy_entries!(J, lop.J)
  J
end

###################################
# NonlinearSolver / StageOperator #
###################################
# Default behaviour from Gridap.Algebra.

#########################################
# NonlinearSolver / LinearStageOperator #
#########################################
# Skip numerical setup update if possible. Since we cannot dispatch on the
# numerical setup to prevent it from updating itself when the matrix is
# constant, we have to overwrite the `NonlinearSolver` interface.

function Algebra._update_nlsolve_cache!(
  cache::NLSolversCache,
  x0::AbstractVector, lop::LinearStageOperator
)
  f!(r, x) = residual!(r, lop, x)
  j!(j, x) = jacobian!(j, lop, x)
  fj!(r, j, x) = residual_and_jacobian!(r, j, lop, x)
  f0, j0 = cache.f0, cache.j0
  residual_and_jacobian!(f0, j0, lop, x0)
  df = OnceDifferentiable(f!, j!, fj!, x0, f0, j0)

  ns = cache.ns
  if !lop.reuse
    numerical_setup!(ns, j0)
  end

  NLSolversCache(f0, j0, df, ns, nothing)
end

function Algebra._nlsolve_with_updated_cache!(
  x::AbstractVector,
  nls::NLSolver, lop::LinearStageOperator,
  cache::NLSolversCache
)
  # After checking NLsolve.jl, the linsolve argument is only passed to Newton
  # and it is only called on the jacobian (j!), so we can save from updating
  # the numerical setup
  ns = cache.ns
  function linsolve!(x, A, b)
    if !lop.reuse
      numerical_setup!(ns, A)
    end
    solve!(x, ns, b)
  end

  result = nlsolve(cache.df, x; linsolve=linsolve!, nls.kwargs...)
  cache.result = result
  copy_entries!(x, result.zero)
end

# IMPORTANT: because `NewtonRaphsonSolver` calls `numerical_setup!` internally,
# we would need to rewrite the functions `solve!` and `_solve_nr!` entirely
# for `LinearStageOperator` in order to skip numerical setup updates when the
# matrix is constant. To be on the safe side, and since `NewtonRaphsonSolver`
# is not exported anyway, we just prevent the user from using it at as a
# nonlinear solver for `LinearStageOperator`.

const nr_on_lop_msg = """
You are trying to use `NewtonRaphsonSolver` to solve a `LinearStageOperator`.
Since this is not optimised (yet), it is forbidden for now. Consider using a
nonlinear solver coming from `NLSolvers`, e.g.
```
  ls = LUSolver()
  nls = NLSolver(ls, show_trace=true, method=:newton, iterations=10)
```
"""

function Algebra.solve!(
  x::AbstractVector,
  nls::NewtonRaphsonSolver, lop::LinearStageOperator,
  cache::Nothing
)
  @unreachable nr_on_lop_msg
end

function Algebra.solve!(
  x::AbstractVector,
  nls::NewtonRaphsonSolver, lop::LinearStageOperator,
  cache
)
  @unreachable nr_on_lop_msg
end

################################
# LinearSolver / StageOperator #
################################
# Forbid solving `StageOperator`s with `LinearSolver`s. For now it is already
# forbidden to solve a generic `FEOperator` with a `LinearSolver`, but it is
# still possible to solve a `NonlinearOperator` with a `LinearSolver`. The
# following should be replicated in Gridap.Algebra for `NonlinearOperator`s at
# some point.
const ls_on_nlop = """
Cannot solve a generic `StageOperator` with a `LinearSolver`.
"""

function Algebra.solve!(
  x::AbstractVector,
  ls::LinearSolver, nlop::StageOperator,
  cache::Nothing
)
  @unreachable ls_on_nlop
end

function Algebra.solve!(
  x::AbstractVector,
  ls::LinearSolver, nlop::StageOperator,
  cache
)
  @unreachable ls_on_nlop
end

######################################
# LinearSolver / LinearStageOperator #
######################################
function Algebra.solve!(
  x::AbstractVector,
  ls::LinearSolver, lop::LinearStageOperator,
  ns::Nothing
)
  J = lop.J
  ss = symbolic_setup(ls, J)
  ns = numerical_setup(ss, J)

  r = lop.r
  rmul!(r, -1)

  solve!(x, ns, r)
  ns
end

function Algebra.solve!(
  x::AbstractVector,
  ls::LinearSolver, lop::LinearStageOperator,
  ns
)
  if !lop.reuse
    J = lop.J
    numerical_setup!(ns, J)
  end

  r = lop.r
  rmul!(r, -1)

  solve!(x, ns, r)
  ns
end
