#######################
# TransientFESolution #
#######################
"""
    struct TransientFESolution <: GridapType end

Wrapper around a `TransientFEOperator` and `ODESolver` that represents the
solution at a set of time steps. It is an iterator that computes the solution
at each time step in a lazy fashion when accessing the solution.

# Mandatory
- [`Base.iterate(fesltn)`](@ref)
- [`Base.iterate(fesltn, state)`](@ref)
"""
struct TransientFESolution <: GridapType
  odesltn::ODESolution
  trial
end

# Constructors
function TransientFESolution(
  odeslvr::ODESolver, feop::TransientFEOperator,
  uh0::Tuple{Vararg{CellField}}, t0::Real, tF::Real
)
  u0 = ()
  for uhi in uh0
    u0 = (u0..., get_free_dof_values(uhi))
  end
  odeop = get_algebraic_operator(feop)
  odesltn = solve(odeslvr, odeop, u0, t0, tF)
  trial = get_trial(feop)
  TransientFESolution(odesltn, trial)
end

function TransientFESolution(
  odeslvr::ODESolver, feop::TransientFEOperator,
  uh0::CellField, t0::Real, tF::Real
)
  TransientFESolution(odeslvr, feop, (uh0,), t0, tF)
end

"""
    solve(
      odeslvr::ODESolver, feop::TransientFEOperator,
      uh0, t0::Real, tF::Real
    ) -> TransientFESolution

Create a `TransientFESolution` wrapper around the `TransientFEOperator` and
`ODESolver`, starting with state `us0` at time `t0`, to be evolved until `tF`.
"""
function Algebra.solve(
  odeslvr::ODESolver, feop::TransientFEOperator,
  uh0, t0::Real, tF::Real
)
  TransientFESolution(odeslvr, feop, uh0, t0, tF)
end

# @fverdugo this is a general implementation of iterate for TransientFESolution
# We could also implement another one for the very common case that the
# underlying odeop is a ODEOpFromFEOp object
function Base.iterate(fesltn::TransientFESolution)
  odesltn_next = iterate(fesltn.odesltn)
  if isnothing(odesltn_next)
    return nothing
  end

  (uF, tF), odesltn_state = odesltn_next

  Uh = allocate_space(fesltn.trial)
  Uh = evaluate!(Uh, fesltn.trial, tF)
  uh = FEFunction(Uh, uF)

  state = (Uh, odesltn_state)
  (uh, tF), state
end

function Base.iterate(fesltn::TransientFESolution, state)
  Uh, odesltn_state = state

  odesltn_next = iterate(fesltn.odesltn, odesltn_state)
  if isnothing(odesltn_next)
    return nothing
  end

  (uF, tF), odesltn_state = odesltn_next

  Uh = evaluate!(Uh, fesltn.trial, tF)
  uh = FEFunction(Uh, uF)

  state = (Uh, odesltn_state)
  (uh, tF), state
end

Base.IteratorSize(::Type{TransientFESolution}) = Base.SizeUnknown()

########
# Test #
########
"""
    test_transient_fe_solution(fesltn::TransientFESolution) -> Bool

Test the interface of `TransientFESolution` specializations.
"""
function test_transient_fe_solution(fesltn::TransientFESolution)
  for (uh_n, t_n) in fesltn
    @test t_n isa Real
    @test uh_n isa FEFunction
  end
  true
end

"""
    test_transient_fe_solver(
      odeslvr::ODESolver, feop::TransientFEOperator,
      uh0, t0, tF
    ) -> Bool

Test the interface of `ODESolver` specializations on `TransientFEOperator`s.
"""
function test_transient_fe_solver(
  odeslvr::ODESolver, feop::TransientFEOperator,
  uh0, t0, tF
)
  fesltn = solve(odeslvr, feop, uh0, t0, tF)
  test_transient_fe_solution(fesltn)
end
