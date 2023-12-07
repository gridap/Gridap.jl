# TODO There is probably more work to do here to implement necessary interface
# functions of FESpace for AbstractTransientTrialFESpace and especially for
# TransientMultiFieldFESpace to implement the MultiFieldFESpace interface

#################################
# AbstractTransientTrialFESpace #
#################################
"""
    abstract type AbstractTransientTrialFESpace <: FESpace end

Transient version of `TrialFESpace`: the Dirichlet boundary conditions are
allowed to be time-dependent.

# Mandatory (in addition to the `FESpace` interface)
- [`allocate_space(space)`](@ref)
- [`evaluate(space, t::Nothing)`](@ref)
- [`evaluate!(space, t)`](@ref)
- [`time_derivative(space)`](@ref)

# Optional
- [`evaluate(space, t::Real)`](@ref)
"""
abstract type AbstractTransientTrialFESpace <: FESpace end

"""
    allocate_space(space::AbstractTransientTrialFESpace) -> FESpace

Allocate a transient space, intended to be updated at every time step.
"""
function allocate_space(space::AbstractTransientTrialFESpace)
  @abstractmethod
end

allocate_space(space::FESpace) = space

"""
    evaluate(space::AbstractTransientTrialFESpace, t::Real) -> FESpace

Allocate a transient space and evaluate the Dirichlet values at time `t`.
"""
function Arrays.evaluate(space::AbstractTransientTrialFESpace, t::Real)
  transient_space = allocate_space(space)
  evaluate!(transient_space, space, t)
  transient_space
end

Arrays.evaluate(space::FESpace, t::Real) = space

"""
    evaluate(space::AbstractTransientTrialFESpace, t::Nothing) -> FESpace

Evaluating at `nothing` means that the Dirichlet values are not important.
"""
function Arrays.evaluate(space::AbstractTransientTrialFESpace, t::Nothing)
  @abstractmethod
end

Arrays.evaluate(space::FESpace, t::Nothing) = space

"""
    evaluate!(
      transient_space::FESpace,
      space::AbstractTransientTrialFESpace, t::Real
    ) -> FESpace

Replace the Dirichlet values of the space by those at time `t`.
"""
function Arrays.evaluate!(
  transient_space::FESpace,
  space::AbstractTransientTrialFESpace, t::Real
)
  @abstractmethod
end

Arrays.evaluate!(transient_space::FESpace, space::FESpace, t::Real) = space

"""
    (space::AbstractTransientTrialFESpace)(t) -> FESpace

Alias for [`evaluate(space, t)`](@ref).
"""
(space::AbstractTransientTrialFESpace)(t) = evaluate(space, t)

# TODO why is this needed?
@static if VERSION >= v"1.3"
  (space::FESpace)(t) = evaluate(space, t)
end

(space::TrialFESpace)(t) = evaluate(space, t)
(space::ZeroMeanFESpace)(t) = evaluate(space, t)

"""
    time_derivative(space::AbstractTransientTrialFESpace) -> FESpace

First-order time derivative of the Dirichlet functions.
"""
function time_derivative(space::AbstractTransientTrialFESpace)
  @abstractmethod
end

# Specialisation for `SingleFieldFESpace`
function time_derivative(space::SingleFieldFESpace)
  HomogeneousTrialFESpace(space)
end

function time_derivative(space::SingleFieldFESpace, ::Val{0})
  space
end

function time_derivative(space::SingleFieldFESpace, ::Val)
  time_derivative(space)
end

# Specialisation for `MultiFieldFESpace`
function time_derivative(space::MultiFieldFESpace)
  MultiFieldFESpace(time_derivative.(space.spaces))
end

# FESpace interface
function FESpaces.SparseMatrixAssembler(
  trial::AbstractTransientTrialFESpace,
  test::FESpace
)
  SparseMatrixAssembler(evaluate(trial, nothing), test)
end

#########################
# TransientTrialFESpace #
#########################
"""
    struct TransientTrialFESpace <: AbstractTransientTrialFESpace end

Single-field `FESpace` with transient Dirichlet data.
"""
struct TransientTrialFESpace{U,U0} <: AbstractTransientTrialFESpace
  space::U
  homogeneous_space::U0
  transient_dirichlet::Union{Function,Tuple{Vararg{Function}}}

  function TransientTrialFESpace(
    space::FESpace, transient_dirichlet::Union{Function,Tuple{Vararg{Function}}}
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

# AbstractTransientTrialFESpace interface
function allocate_space(U::TransientTrialFESpace)
  HomogeneousTrialFESpace(U.space)
end

Arrays.evaluate(U::TransientTrialFESpace, t::Nothing) = U.homogeneous_space

function Arrays.evaluate!(Ut::FESpace, U::TransientTrialFESpace, t::Real)
  if U.transient_dirichlet isa AbstractVector
    dirichlets_at_t = map(o -> o(t), U.transient_dirichlet)
  else
    dirichlets_at_t = U.transient_dirichlet(t)
  end
  TrialFESpace!(Ut, dirichlets_at_t)
  Ut
end

function time_derivative(U::TransientTrialFESpace)
  TransientTrialFESpace(U.space, time_derivative.(U.transient_dirichlet))
end

# FESpace interface
FESpaces.zero_free_values(f::TransientTrialFESpace) = zero_free_values(f.space)

FESpaces.get_vector_type(f::TransientTrialFESpace) = get_vector_type(f.space)

FESpaces.get_dof_value_type(f::TransientTrialFESpace) = get_dof_value_type(f.space)

FESpaces.has_constraints(f::TransientTrialFESpace) = has_constraints(f.space)

##############################
# TransientMultiFieldFESpace #
##############################
# (Copy-pasted from MultiField)
"""
    struct TransientMultiFieldFESpace <: AbstractTransientTrialFESpace

Multi-field `FESpace` with transient Dirichlet data.
"""
struct TransientMultiFieldFESpace{
  MS<:MultiFieldStyle,CS<:ConstraintStyle,V
} <: AbstractTransientTrialFESpace
  vector_type::Type{V}
  spaces::Vector
  multi_field_style::MS
  constraint_style::CS

  function TransientMultiFieldFESpace(
    ::Type{V}, spaces::Vector,
    multi_field_style::MultiFieldStyle
  ) where {V}
    @assert length(spaces) > 0

    MS = typeof(multi_field_style)
    if any(map(has_constraints, spaces))
      constraint_style = Constrained()
    else
      constraint_style = UnConstrained()
    end
    CS = typeof(constraint_style)
    new{MS,CS,V}(V, spaces, multi_field_style, constraint_style)
  end
end

# Constructors (copy-pasted from MultiField)
function TransientMultiFieldFESpace(
  spaces::Vector;
  style=ConsecutiveMultiFieldStyle()
)
  Ts = map(get_dof_value_type, spaces)
  T = typeof(*(map(zero, Ts)...))
  if style isa BlockMultiFieldStyle
    style = BlockMultiFieldStyle(style, spaces)
    VT = typeof(mortar(map(zero_free_values, spaces)))
  else
    VT = Vector{T}
  end
  TransientMultiFieldFESpace(VT, spaces, style)
end

function TransientMultiFieldFESpace(
  ::Type{V}, spaces::AbstractVector{<:FESpace}
) where {V}
  TransientMultiFieldFESpace(V, spaces, ConsecutiveMultiFieldStyle())
end

# Specialisations when none of the spaces is transient
function TransientMultiFieldFESpace(
  spaces::Vector{<:SingleFieldFESpace};
  style=ConsecutiveMultiFieldStyle()
)
  MultiFieldFESpace(spaces; style)
end

function TransientMultiFieldFESpace(
  ::Type{V}, spaces::AbstractVector{<:SingleFieldFESpace}
) where {V}
  MultiFieldFESpace(V, spaces, ConsecutiveMultiFieldStyle())
end

# FESpace interface (copy-pasted from MultiField)
function FESpaces.zero_free_values(f::TransientMultiFieldFESpace{<:BlockMultiFieldStyle{NB,SB,P}}) where {NB,SB,P}
  block_ranges = get_block_ranges(NB, SB, P)
  block_num_dofs = map(range -> sum(map(num_free_dofs, f.spaces[range])), block_ranges)
  block_vtypes = map(range -> get_vector_type(first(f.spaces[range])), block_ranges)
  return mortar(map(allocate_vector, block_vtypes, block_num_dofs))
end

FESpaces.get_dof_value_type(f::TransientMultiFieldFESpace{MS,CS,V}) where {MS,CS,V} = eltype(V)

FESpaces.get_vector_type(f::TransientMultiFieldFESpace) = f.vector_type

FESpaces.ConstraintStyle(::Type{TransientMultiFieldFESpace{S,B,V}}) where {S,B,V} = B()
FESpaces.ConstraintStyle(::TransientMultiFieldFESpace) = ConstraintStyle(typeof(f))

# Overwrite SparseMatrixAssembler to use BlockSparseMatrixAssembler instead of
# GenericSparseMatrixAssembler
function FESpaces.SparseMatrixAssembler(
  mat, vec,
  trial::TransientMultiFieldFESpace{MS},
  test::TransientMultiFieldFESpace{MS},
  strategy::AssemblyStrategy=DefaultAssemblyStrategy()
) where {MS<:BlockMultiFieldStyle}
  style = MultiFieldStyle(test)
  mat_builder = SparseMatrixBuilder(mat)
  vec_builder = ArrayBuilder(vec)

  BlockSparseMatrixAssembler(
    style, trial, test,
    mat_builder, vec_builder,
    strategy
  )
end

# MultiField interface (copy-pasted from MultiField)
MultiField.MultiFieldStyle(::Type{TransientMultiFieldFESpace{S,B,V}}) where {S,B,V} = S()
MultiField.MultiFieldStyle(f::TransientMultiFieldFESpace) = MultiFieldStyle(typeof(f))

function MultiField.num_fields(f::TransientMultiFieldFESpace)
  length(f.spaces)
end

Base.iterate(m::TransientMultiFieldFESpace) = iterate(m.spaces)

Base.iterate(m::TransientMultiFieldFESpace, state) = iterate(m.spaces, state)

Base.getindex(m::TransientMultiFieldFESpace, field_id::Integer) = m.spaces[field_id]

Base.length(m::TransientMultiFieldFESpace) = length(m.spaces)

# AbstractTransientTrialFESpace interface
function allocate_space(U::TransientMultiFieldFESpace)
  spaces = [allocate_space(Ui) for Ui in U]
  style = MultiFieldStyle(U)
  MultiFieldFESpace(spaces; style)
end

function Arrays.evaluate(U::TransientMultiFieldFESpace, t::Nothing)
  spaces = [evaluate(Ui, nothing) for Ui in U]
  style = MultiFieldStyle(U)
  MultiFieldFESpace(spaces; style)
end

function Arrays.evaluate!(Ut::FESpace, U::TransientMultiFieldFESpace, t::Real)
  spaces = [evaluate!(Uti, Ui, t) for (Uti, Ui) in zip(Ut, U)]
  style = MultiFieldStyle(U)
  MultiFieldFESpace(spaces; style)
end

function time_derivative(U::TransientMultiFieldFESpace)
  spaces = time_derivative.(U.spaces)
  style = MultiFieldStyle(U)
  TransientMultiFieldFESpace(spaces; style)
end

########
# Test #
########
"""
    test_transient_trial_fe_space(U::AbstractTransientTrialFESpace) -> Bool

Test the interface of `AbstractTransientTrialFESpace` specializations.
"""
function test_transient_trial_fe_space(U::AbstractTransientTrialFESpace)
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
