#########################
# TransientTrialFESpace #
#########################
"""
    struct TransientTrialFESpace <: SingleFieldFESpace end

Transient version of `TrialFESpace`: the Dirichlet boundary conditions are
allowed to be time-dependent.

# Mandatory
- [`allocate_space(space)`](@ref)
- [`evaluate!(space, t)`](@ref)
- [`evaluate(space, t)`](@ref)
- [`time_derivative(space)`](@ref)

# Optional
- [`evaluate(space, t::Real)`](@ref)
"""
struct TransientTrialFESpace{U,U0} <: SingleFieldFESpace
  space::U
  homogeneous_space::U0
  transient_dirichlet::Union{Function,AbstractVector{<:Function}}

  function TransientTrialFESpace(
    space::FESpace, transient_dirichlet::Union{Function,AbstractVector{<:Function}}
  )
    homogeneous_space = HomogeneousTrialFESpace(space)
    U = typeof(space)
    U0 = typeof(homogeneous_space)
    new{U,U0}(space, homogeneous_space, transient_dirichlet)
  end
end

# Constructors
function TransientTrialFESpace(space)
  HomogeneousTrialFESpace(space)
end

"""
    allocate_space(space::TransientTrialFESpace) -> FESpace

Allocate a transient space, intended to be updated at every time step.
"""
function allocate_space(U::TransientTrialFESpace)
  HomogeneousTrialFESpace(U.space)
end

"""
    evaluate!(
      transient_space::FESpace,
      space::TransientTrialFESpace, t::Real
    ) -> FESpace

Replace the Dirichlet values of the space by those at time `t`.
"""
function Arrays.evaluate!(Ut::FESpace, U::TransientTrialFESpace, t::Real)
  if U.transient_dirichlet isa AbstractVector
    dirichlets_at_t = map(o -> o(t), U.transient_dirichlet)
  else
    dirichlets_at_t = U.transient_dirichlet(t)
  end
  TrialFESpace!(Ut, dirichlets_at_t)
  Ut
end

"""
    evaluate(space::TransientTrialFESpace, t::Real) -> FESpace

Allocate a transient space and evaluate the Dirichlet values at time `t`.
"""
function Arrays.evaluate(U::TransientTrialFESpace, t::Real)
  Ut = allocate_space(U)
  evaluate!(Ut, U, t)
  Ut
end

"""
    evaluate(space::TransientTrialFESpace, t::Nothing) -> FESpace

Evaluating at `nothing` means that the Dirichlet values are not important.
"""
function Arrays.evaluate(U::TransientTrialFESpace, t::Nothing)
  U.homogeneous_space
end

"""
    (space::TransientTrialFESpace)(t) -> FESpace

Alias for [`evaluate(space, t)`](@ref).
"""
(space::TransientTrialFESpace)(t) = evaluate(space, t)

"""
    time_derivative(space::TransientTrialFESpace) -> FESpace

First-order time derivative of the Dirichlet functions.
"""
function time_derivative(U::TransientTrialFESpace)
  TransientTrialFESpace(U.space, time_derivative.(U.transient_dirichlet))
end

# FESpace interface
FESpaces.get_free_dof_ids(f::TransientTrialFESpace) = get_free_dof_ids(f.space)
FESpaces.get_vector_type(f::TransientTrialFESpace) = get_vector_type(f.space)
Geometry.get_triangulation(f::TransientTrialFESpace) = get_triangulation(f.space)
FESpaces.get_cell_dof_ids(f::TransientTrialFESpace) = get_cell_dof_ids(f.space)
FESpaces.get_fe_basis(f::TransientTrialFESpace) = get_fe_basis(f.space)
FESpaces.get_fe_dof_basis(f::TransientTrialFESpace) = get_fe_dof_basis(f.space)
FESpaces.ConstraintStyle(::Type{<:TransientTrialFESpace{U}}) where {U} = ConstraintStyle(U)
function FESpaces.get_cell_constraints(f::TransientTrialFESpace, c::Constrained)
  get_cell_constraints(f.space, c)
end
function FESpaces.get_cell_isconstrained(f::TransientTrialFESpace, c::Constrained)
  get_cell_isconstrained(f.space, c)
end

# SingleFieldFESpace interface
FESpaces.get_dirichlet_dof_ids(f::TransientTrialFESpace) = get_dirichlet_dof_ids(f.space)
FESpaces.num_dirichlet_tags(f::TransientTrialFESpace) = num_dirichlet_tags(f.space)
FESpaces.get_dirichlet_dof_tag(f::TransientTrialFESpace) = get_dirichlet_dof_tag(f.space)
function FESpaces.scatter_free_and_dirichlet_values(f::TransientTrialFESpace, free_values, dirichlet_values)
  scatter_free_and_dirichlet_values(f.space, free_values, dirichlet_values)
end
function FESpaces.gather_free_and_dirichlet_values!(free_values, dirichlet_values, f::TransientTrialFESpace, cell_vals)
  gather_free_and_dirichlet_values!(free_values, dirichlet_values, f.space, cell_vals)
end

function FESpaces.get_dirichlet_dof_values(f::TransientTrialFESpace)
  msg = """
  It does not make sense to get the Dirichlet DOF values of a transient FE space. You
  should first evaluate the transient FE space at a point in time and get the Dirichlet
  DOF values from there.
  """
  @unreachable msg
end

# function FESpaces.SparseMatrixAssembler(
#   trial::TransientTrialFESpace,
#   test::FESpace
# )
#   SparseMatrixAssembler(evaluate(trial, nothing), test)
# end

###########
# FESpace #
###########
allocate_space(space::FESpace) = space
Arrays.evaluate!(transient_space::FESpace, space::FESpace, t::Real) = space
Arrays.evaluate(space::FESpace, t::Real) = space
Arrays.evaluate(space::FESpace, t::Nothing) = space

# TODO why is this needed?
@static if VERSION >= v"1.3"
  (space::FESpace)(t) = evaluate(space, t)
end
(space::TrialFESpace)(t) = evaluate(space, t)
(space::ZeroMeanFESpace)(t) = evaluate(space, t)

function time_derivative(space::SingleFieldFESpace)
  HomogeneousTrialFESpace(space)
end

#####################
# MultiFieldFESpace #
#####################
# This is only for backward compatibility, we could remove it
const TransientMultiFieldFESpace = MultiFieldFESpace

function has_transient(U::MultiFieldFESpace)
  any(space -> space isa TransientTrialFESpace, U.spaces)
end

function allocate_space(U::MultiFieldFESpace)
  if !has_transient(U)
    return U
  end
  spaces = map(allocate_space, U)
  style = MultiFieldStyle(U)
  MultiFieldFESpace(spaces; style)
end

function Arrays.evaluate!(Ut::MultiFieldFESpace, U::MultiFieldFESpace, t::Real)
  if !has_transient(U)
    return Ut
  end
  for (Uti, Ui) in zip(Ut, U)
    evaluate!(Uti, Ui, t)
  end
  Ut
end

function Arrays.evaluate(U::MultiFieldFESpace, t::Real)
  if !has_transient(U)
    return U
  end
  Ut = allocate_space(U)
  evaluate!(Ut, U, t)
end

function Arrays.evaluate(U::MultiFieldFESpace, t::Nothing)
  if !has_transient(U)
    return U
  end
  spaces = map(space -> evaluate(space, t), U.spaces)
  style = MultiFieldStyle(U)
  MultiFieldFESpace(spaces; style)
end

function time_derivative(U::MultiFieldFESpace)
  spaces = map(time_derivative, U.spaces)
  style = MultiFieldStyle(U)
  MultiFieldFESpace(spaces; style)
end

########
# Test #
########
"""
    test_tfe_space(U::FESpace) -> Bool

Test the transient interface of `FESpace` specializations.
"""
function test_tfe_space(U::FESpace)
  UX = evaluate(U, nothing)
  @test UX isa FESpace

  t = 0.0

  U0 = allocate_space(U)
  U0 = evaluate!(U0, U, t)
  @test U0 isa FESpace

  U0 = evaluate(U, t)
  @test U0 isa FESpace

  U0 = U(t)
  @test U0 isa FESpace

  Ut = ∂t(U)
  Ut0 = Ut(t)
  @test Ut0 isa FESpace

  true
end
