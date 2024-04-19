#######################
# TransientFESolution #
#######################
"""
    abstract type TransientFESolution <: GridapType end

Wrapper around a `TransientFEOperator` and `ODESolver` that represents the
solution at a set of time steps. It is an iterator that computes the solution
at each time step in a lazy fashion when accessing the solution.

# Mandatory
- [`iterate(tfesltn)`](@ref)
- [`iterate(tfesltn, state)`](@ref)
"""
abstract type TransientFESolution <: GridapType end

"""
    Base.iterate(tfesltn::TransientFESolution) -> ((Real, FEFunction), StateType)

Allocate a cache and perform one step of the `ODEOperator` with the `ODESolver`
attached to the `TransientFESolution`.
"""
function Base.iterate(tfesltn::TransientFESolution)
  @abstractmethod
end

"""
    Base.iterate(tfesltn::TransientFESolution) -> ((Real, FEFunction), StateType)

Perform one step of the `ODEOperator` with the `ODESolver` attached to the
`TransientFESolution`.
"""
function Base.iterate(tfesltn::TransientFESolution, state)
  @abstractmethod
end

Base.IteratorSize(::Type{<:TransientFESolution}) = Base.SizeUnknown()

##############################
# GenericTransientFESolution #
##############################
"""
    struct GenericTransientFESolution <: TransientFESolution end

Generic wrapper for the evolution of an `TransientFEOperator` with an
`ODESolver`.
"""

struct GenericTransientFESolution <: TransientFESolution
  odesltn::ODESolution
  trial
end

# Constructors
function GenericTransientFESolution(
  odeslvr::ODESolver, tfeop::TransientFEOperator,
  t0::Real, tF::Real, uhs0::Tuple{Vararg{CellField}}
)
  odeop = get_algebraic_operator(tfeop)
  us0 = get_free_dof_values.(uhs0)
  odesltn = solve(odeslvr, odeop, t0, tF, us0)
  trial = get_trial(tfeop)
  GenericTransientFESolution(odesltn, trial)
end

function GenericTransientFESolution(
  odeslvr::ODESolver, tfeop::TransientFEOperator,
  t0::Real, tF::Real, uh0::CellField,
)
  uhs0 = (uh0,)
  GenericTransientFESolution(odeslvr, tfeop, t0, tF, uhs0)
end

function Base.iterate(tfesltn::GenericTransientFESolution)
  ode_it = iterate(tfesltn.odesltn)
  if isnothing(ode_it)
    return nothing
  end

  ode_it_data, ode_it_state = ode_it
  tF, uF = ode_it_data

  Uh = allocate_space(tfesltn.trial)
  Uh = evaluate!(Uh, tfesltn.trial, tF)
  uhF = FEFunction(Uh, uF)

  tfe_it_data = (tF, uhF)
  tfe_it_state = (Uh, ode_it_state)
  (tfe_it_data, tfe_it_state)
end

function Base.iterate(tfesltn::GenericTransientFESolution, state)
  Uh, ode_it_state = state

  ode_it = iterate(tfesltn.odesltn, ode_it_state)
  if isnothing(ode_it)
    return nothing
  end

  ode_it_data, ode_it_state = ode_it
  tF, uF = ode_it_data

  Uh = evaluate!(Uh, tfesltn.trial, tF)
  uhF = FEFunction(Uh, uF)

  tfe_it_data = (tF, uhF)
  tfe_it_state = (Uh, ode_it_state)
  (tfe_it_data, tfe_it_state)
end

##############################
# Default behaviour of solve #
##############################
"""
    solve(
      odeslvr::ODESolver, tfeop::TransientFEOperator,
      t0::Real, tF::Real, uhs0
    ) -> TransientFESolution

Create a `TransientFESolution` wrapper around the `TransientFEOperator` and
`ODESolver`, starting at time `t0` with state `us0`, to be evolved until `tF`.
"""
function Algebra.solve(
  odeslvr::ODESolver, tfeop::TransientFEOperator,
  t0::Real, tF::Real, uhs0::Tuple{Vararg{CellField}}
)
  GenericTransientFESolution(odeslvr, tfeop, t0, tF, uhs0)
end

function Algebra.solve(
  odeslvr::ODESolver, tfeop::TransientFEOperator,
  t0::Real, tF::Real, uh0::CellField
)
  uhs0 = (uh0,)
  solve(odeslvr, tfeop, t0, tF, uhs0)
end

########
# Test #
########
"""
    test_tfe_solution(tfesltn::TransientFESolution) -> Bool

Test the interface of `TransientFESolution` specializations.
"""
function test_tfe_solution(tfesltn::TransientFESolution)
  for (t_n, uh_n) in tfesltn
    @test t_n isa Real
    @test uh_n isa FEFunction
  end

  true
end

"""
    test_tfe_solver(
      odeslvr::ODESolver, tfeop::TransientFEOperator,
      t0::Real, tF::Real, uhs0
    ) -> Bool

Test the interface of `ODESolver` specializations on `TransientFEOperator`s.
"""
function test_tfe_solver(
  odeslvr::ODESolver, tfeop::TransientFEOperator,
  t0::Real, tF::Real, uhs0::Tuple{Vararg{AbstractVector}}
)
  tfesltn = solve(odeslvr, tfeop, t0, tF, uhs0)
  test_tfe_solution(tfesltn)
end
