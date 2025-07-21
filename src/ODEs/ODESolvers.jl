#############
# ODESolver #
#############
"""
    abstract type ODESolver <: GridapType end

An `ODESolver` is a map that update state vectors. These state vectors are
created at the first iteration from the initial conditions, and are then
converted back into the evaluation of the solution at the current time step.

In the simplest case, the state vectors correspond to the first `N-1` time
derivatives of `u` at time `t_n`, where `N` is the order of the `ODEOperator`,
but some solvers rely on other state variables (values at previous times,
 higher-order derivatives...).

# Mandatory
- [`allocate_odecache(odeslvr, odeop, t0, us0)`](@ref)
- [`ode_march!(stateF, odeslvr, odeop, t0, state0, odecache)`](@ref)

# Optional
- [`ode_start(odeslvr, odeop, t0, us0, odecache)`](@ref)
- [`ode_finish!(uF, odeslvr, odeop, t0, tF, stateF, odecache)`](@ref)
"""
abstract type ODESolver <: GridapType end

"""
    allocate_odecache(
      odeslvr::ODESolver, odeop::ODEOperator,
      t0::Real, us0::Tuple{Vararg{AbstractVector}}
    ) -> CacheType

Allocate the cache of the `ODESolver` applied to the `ODEOperator`.
"""
function allocate_odecache(
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, us0::Tuple{Vararg{AbstractVector}}
)
  @abstractmethod
end

"""
    ode_start(
      odeslvr::ODESolver, odeop::ODEOperator,
      t0::Real, us0::Tuple{Vararg{AbstractVector}},
      odecache
    ) -> (Tuple{Vararg{AbstractVector}}, CacheType)

Convert the initial conditions into state vectors.
"""
function ode_start(
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, us0::Tuple{Vararg{AbstractVector}},
  odecache
)
  state0 = copy.(us0)
  (state0, odecache)
end

"""
    ode_march!(
      stateF::Tuple{Vararg{AbstractVector}},
      odeslvr::ODESolver, odeop::ODEOperator,
      t0::Real, state0::Tuple{Vararg{AbstractVector}},
      odecache
    ) -> (Real, Tuple{Vararg{AbstractVector}}, CacheType)

March the state vector for one time step.
"""
function ode_march!(
  stateF::Tuple{Vararg{AbstractVector}},
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, state0::Tuple{Vararg{AbstractVector}},
  odecache
)
  @abstractmethod
end

"""
    ode_finish!(
      uF::AbstractVector,
      odeslvr::ODESolver, odeop::ODEOperator,
      t0::Real, tF, stateF::Tuple{Vararg{AbstractVector}},
      odecache
    ) -> (AbstractVector, CacheType)

Convert the state vectors into the evaluation of the solution of the ODE at the
current time.
"""
function ode_finish!(
  uF::AbstractVector,
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, tF, stateF::Tuple{Vararg{AbstractVector}},
  odecache
)
  copy!(uF, first(stateF))
  (uF, odecache)
end

########
# Test #
########
"""
    test_ode_solver(
      odeslvr::ODESolver, odeop::ODEOperator,
      t0::Real, us0::Tuple{Vararg{AbstractVector}}
    ) -> Bool

Test the interface of `ODESolver` specializations.
"""
function test_ode_solver(
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, us0::Tuple{Vararg{AbstractVector}}
)
  odecache = allocate_odecache(odeslvr, odeop, t0, us0)

  # Starting procedure
  state0, odecache = ode_start(
    odeslvr, odeop,
    t0, us0,
    odecache
  )
  @test state0 isa Tuple{Vararg{AbstractVector}}

  # Marching procedure
  stateF = copy.(state0)
  tF, stateF, odecache = ode_march!(
    stateF,
    odeslvr, odeop,
    t0, state0,
    odecache
  )
  @test tF isa Real
  @test stateF isa Tuple{Vararg{AbstractVector}}

  # Finishing procedure
  uF = copy(first(us0))
  uF, odecache = ode_finish!(
    uF,
    odeslvr, odeop,
    t0, tF, stateF,
    odecache
  )
  @test uF isa AbstractVector

  true
end

##################
# Import solvers #
##################
# First-order
include("ODESolvers/ForwardEuler.jl")

include("ODESolvers/ThetaMethod.jl")

include("ODESolvers/GeneralizedAlpha1.jl")

include("ODESolvers/Tableaus.jl")

include("ODESolvers/RungeKuttaEX.jl")

include("ODESolvers/RungeKuttaDIM.jl")

include("ODESolvers/RungeKuttaIMEX.jl")

# Second-order
include("ODESolvers/GeneralizedAlpha2.jl")

#########
# Utils #
#########
function _setindex_all!(a::CompressedArray, v, i::Integer)
  # This is a straightforward implementation of setindex! for `CompressedArray`
  # when we want to update the value associated to all pointers currently
  # pointing to the same value
  idx = a.ptrs[i]
  a.values[idx] = v
  a
end

"""
    RungeKutta(sysslvr_nl::NonlinearSolver, sysslvr_l::NonlinearSolver, dt::Real, tableau::AbstractTableau)
    RungeKutta(sysslvr_nl::NonlinearSolver, dt::Real, tableau)
    RungeKutta(sysslvr_nl::NonlinearSolver, sysslvr_l::NonlinearSolver, dt::Real, name::Symbol)

The second constructor uses `sysslvr_nl` for the `sysslvr_l` argument.
"""
function RungeKutta(
  sysslvr_nl::NonlinearSolver, sysslvr_l::NonlinearSolver,
  dt::Real, tableau::AbstractTableau
)
  type = TableauType(tableau)
  if type == ExplicitTableau
    EXRungeKutta(sysslvr_nl, dt, tableau)
  elseif type == DiagonallyImplicitTableau
    DIMRungeKutta(sysslvr_nl, sysslvr_l, dt, tableau)
  elseif type == ImplicitExplicitTableau
    IMEXRungeKutta(sysslvr_nl, sysslvr_l, dt, tableau)
    # elseif type == FullyImplicitTableau
    #   FIMRungeKutta(sysslvr_nl, sysslvr_l, dt, tableau)
  end
end

function RungeKutta(
  sysslvr_nl::NonlinearSolver, sysslvr_l::NonlinearSolver,
  dt::Real, name::Symbol
)
  RungeKutta(sysslvr_nl, sysslvr_l, dt, ButcherTableau(name))
end

function RungeKutta(sysslvr_nl::NonlinearSolver, dt::Real, tableau)
  RungeKutta(sysslvr_nl, sysslvr_nl, dt, tableau)
end
