#######################
# TransientFESolution #
#######################
"""
    abstract type TransientFESolution <: GridapType end

Wrapper around a `TransientFEOperator` and `ODESolver` that represents the
solution at a set of time steps. It is an iterator that computes the solution
at each time step in a lazy fashion when accessing the solution.

# Mandatory
- [`Base.iterate(fesltn)`](@ref)
- [`Base.iterate(fesltn, state)`](@ref)
"""
abstract type TransientFESolution <: GridapType end

"""
    Base.iterate(fesltn::TransientFESolution) -> ((Real, FEFunction), StateType)

Allocate a cache and perform one step of the `ODEOperator` with the `ODESolver`
attached to the `TransientFESolution`.
"""
function Base.iterate(odesltn::TransientFESolution)
  @abstractmethod
end

"""
    Base.iterate(fesltn::TransientFESolution) -> ((Real, FEFunction), StateType)

Perform one step of the `ODEOperator` with the `ODESolver` attached to the
`TransientFESolution`.
"""
function Base.iterate(odesltn::TransientFESolution, state)
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
  odeslvr::ODESolver, feop::TransientFEOperator,
  t0::Real, tF::Real, uhs0::Tuple{Vararg{CellField}}
)
  us0 = ()
  for uhi0 in uhs0
    ui0 = get_free_dof_values(uhi0)
    us0 = (us0..., ui0)
  end
  odeop = get_algebraic_operator(feop)
  odesltn = solve(odeslvr, odeop, t0, tF, us0)
  trial = get_trial(feop)
  GenericTransientFESolution(odesltn, trial)
end

function GenericTransientFESolution(
  odeslvr::ODESolver, feop::TransientFEOperator,
  t0::Real, tF::Real, uh0::CellField,
)
  GenericTransientFESolution(odeslvr, feop, t0, tF, (uh0,))
end

function Base.iterate(fesltn::GenericTransientFESolution)
  odesltn_next = iterate(fesltn.odesltn)
  if isnothing(odesltn_next)
    return nothing
  end

  (tF, uF), odesltn_state = odesltn_next

  Uh = allocate_space(fesltn.trial)
  Uh = evaluate!(Uh, fesltn.trial, tF)
  uhF = FEFunction(Uh, uF)

  state = (Uh, odesltn_state)
  (tF, uhF), state
end

function Base.iterate(fesltn::GenericTransientFESolution, state)
  Uh, odesltn_state = state

  odesltn_next = iterate(fesltn.odesltn, odesltn_state)
  if isnothing(odesltn_next)
    return nothing
  end

  (tF, uF), odesltn_state = odesltn_next

  Uh = evaluate!(Uh, fesltn.trial, tF)
  uhF = FEFunction(Uh, uF)

  state = (Uh, odesltn_state)
  (tF, uhF), state
end

Base.IteratorSize(::Type{GenericTransientFESolution}) = Base.SizeUnknown()

##############################
# Default behaviour of solve #
##############################
"""
    solve(
      odeslvr::ODESolver, feop::TransientFEOperator,
      t0::Real, tF::Real, uhs0
    ) -> TransientFESolution

Create a `TransientFESolution` wrapper around the `TransientFEOperator` and
`ODESolver`, starting at time `t0` with state `us0`, to be evolved until `tF`.
"""
function Algebra.solve(
  odeslvr::ODESolver, feop::TransientFEOperator,
  t0::Real, tF::Real, uhs0::Tuple{Vararg{CellField}}
)
  GenericTransientFESolution(odeslvr, feop, t0, tF, uhs0)
end

function Algebra.solve(
  odeslvr::ODESolver, feop::TransientFEOperator,
  t0::Real, tF::Real, uh0::CellField
)
  GenericTransientFESolution(odeslvr, feop, t0, tF, (uh0,))
end

########
# Test #
########
"""
    test_transient_fe_solution(fesltn::TransientFESolution) -> Bool

Test the interface of `TransientFESolution` specializations.
"""
function test_transient_fe_solution(fesltn::TransientFESolution)
  for (t_n, uh_n) in fesltn
    @test t_n isa Real
    @test uh_n isa FEFunction
  end

  true
end

"""
    test_transient_fe_solver(
      odeslvr::ODESolver, feop::TransientFEOperator,
      t0::Real, tF::Real, uhs0
    ) -> Bool

Test the interface of `ODESolver` specializations on `TransientFEOperator`s.
"""
function test_transient_fe_solver(
  odeslvr::ODESolver, feop::TransientFEOperator,
  t0::Real, tF::Real, uhs0::Tuple{Vararg{AbstractVector}}
)
  fesltn = solve(odeslvr, feop, t0, tF, uhs0)
  test_transient_fe_solution(fesltn)
end
