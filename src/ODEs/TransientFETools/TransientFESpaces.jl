#########################
# TransientTrialFESpace #
#########################
# QUESTION Should TransientTrialFESpace be a subtype of FESpace? Same question
# for TransientMultiFieldTrialFESpace. Alternatively, could create a
# TransientFESpace subtype of FESpace.
"""
A single field FE space with transient Dirichlet data (see Multifield below).
"""
struct TransientTrialFESpace{A,B}
  space::A
  dirichlet_t::Union{Function,AbstractVector{<:Function}}
  Ud0::B

  function TransientTrialFESpace(space::A,
    dirichlet_t::Union{Function,AbstractVector{<:Function}}) where {A}
    Ud0 = HomogeneousTrialFESpace(space)
    B = typeof(Ud0)
    new{A,B}(space, dirichlet_t, Ud0)
  end
end

function TransientTrialFESpace(space::A) where {A}
  HomogeneousTrialFESpace(space)
end

# Update of the Dirichlet values
"""
Allocate the space to be used as first argument in evaluate!
"""
function allocate_trial_space(U::TransientTrialFESpace)
  HomogeneousTrialFESpace(U.space)
end

allocate_trial_space(U::FESpace) = U

"""
Allocate a trial space and evaluate the Dirichlet values
"""
function evaluate(U::TransientTrialFESpace, t::Real)
  Ut = allocate_trial_space(U)
  evaluate!(Ut, U, t)
  Ut
end

evaluate(U::FESpace, t::Real) = U

"""
Update the Dirichlet values in place
"""
function evaluate!(Ut::T, U::TransientTrialFESpace, t::Real) where {T}
  if U.dirichlet_t isa AbstractVector
    objects_at_t = map(o -> o(t), U.dirichlet_t)
  else
    objects_at_t = U.dirichlet_t(t)
  end
  TrialFESpace!(Ut, objects_at_t)
  Ut
end

evaluate!(Ut::FESpace, U::FESpace, t::Real) = U

"""
We can evaluate at `nothing` when we do not care about the Dirichlet values
"""
evaluate(U::TransientTrialFESpace, t::Nothing) = U.Ud0
evaluate(U::TrialFESpace, t::Nothing) = U
evaluate(U::FESpace, t::Nothing) = U

"""
Functor-like evaluation. It allocates Dirichlet values in general
"""
(U::TransientTrialFESpace)(t) = evaluate(U, t)
(U::TrialFESpace)(t) = U
(U::ZeroMeanFESpace)(t) = U
# (U::Union{TrialFESpace,ZeroMeanFESpace})(t) = U

@static if VERSION >= v"1.3"
  (U::FESpace)(t) = U
end

# Time derivatives
"""
First time derivative of the Dirichlet functions
"""
∂t(U::TransientTrialFESpace) = TransientTrialFESpace(U.space, ∂t.(U.dirichlet_t))
∂t(U::SingleFieldFESpace) = HomogeneousTrialFESpace(U)
∂t(U::MultiFieldFESpace) = MultiFieldFESpace(∂t.(U.spaces))
∂t(t::T) where {T<:Number} = zero(T)

"""
Second time derivative of the Dirichlet functions
"""
∂tt(U::TransientTrialFESpace) = TransientTrialFESpace(U.space, ∂tt.(U.dirichlet_t))
∂tt(U::SingleFieldFESpace) = HomogeneousTrialFESpace(U)
∂tt(U::MultiFieldFESpace) = MultiFieldFESpace(∂tt.(U.spaces))
∂tt(t::T) where {T<:Number} = zero(T)

# FESpace interface
zero_free_values(f::TransientTrialFESpace) = zero_free_values(f.space)

get_vector_type(f::TransientTrialFESpace) = get_vector_type(f.space)

get_dof_value_type(f::TransientTrialFESpace) = get_dof_value_type(f.space)

has_constraints(f::TransientTrialFESpace) = has_constraints(f.space)

###################################
# TransientMultiFieldTrialFESpace #
###################################
"""
A multi field FE space with transient Dirichlet data (see Multifield below).
"""
struct TransientMultiFieldTrialFESpace{MS<:MultiFieldStyle,CS<:ConstraintStyle,V}
  vector_type::Type{V}
  spaces::AbstractVector
  multi_field_style::MS
  constraint_style::CS

  function TransientMultiFieldTrialFESpace(::Type{V}, spaces::AbstractVector,
    multi_field_style::MultiFieldStyle) where {V}
    @assert length(spaces) > 0

    if any(map(has_constraints, spaces))
      constraint_style = Constrained()
    else
      constraint_style = UnConstrained()
    end

    MS = typeof(multi_field_style)
    CS = typeof(constraint_style)
    new{MS,CS,V}(V, spaces, multi_field_style, constraint_style)
  end
end

# Default constructors
# QUESTION why TransientMultiFieldFESpace instead of
# TransientMultiFieldTrialFESpace?
function TransientMultiFieldFESpace(spaces::AbstractVector;
  style=ConsecutiveMultiFieldStyle())
  Ts = map(get_dof_value_type, spaces)
  T = promote_type(Ts...)

  if style isa BlockMultiFieldStyle
    style = BlockMultiFieldStyle(style, spaces)
    VT = typeof(mortar(map(zero_free_values, spaces)))
  else
    VT = Vector{T}
  end

  TransientMultiFieldTrialFESpace(VT, spaces, style)
end

function TransientMultiFieldFESpace(
  spaces::AbstractVector{<:SingleFieldFESpace};
  style=ConsecutiveMultiFieldStyle())
  MultiFieldFESpace(spaces; style)
end

function TransientMultiFieldFESpace(::Type{V}, spaces::AbstractVector) where {V}
  TransientMultiFieldTrialFESpace(V, spaces, ConsecutiveMultiFieldStyle())
end

function TransientMultiFieldFESpace(::Type{V},
  spaces::AbstractVector{<:SingleFieldFESpace}) where {V}
  MultiFieldFESpace(V, spaces, ConsecutiveMultiFieldStyle())
end

# Base interface
iterate(m::TransientMultiFieldTrialFESpace) = iterate(m.spaces)
iterate(m::TransientMultiFieldTrialFESpace, state) = iterate(m.spaces, state)
getindex(m::TransientMultiFieldTrialFESpace, idx::Integer) = m.spaces[idx]
length(m::TransientMultiFieldTrialFESpace) = length(m.spaces)

# TransientTrialFESpace interface
function allocate_trial_space(U::TransientMultiFieldTrialFESpace)
  spaces = allocate_trial_space.(U.spaces)
  style = MultiFieldStyle(U)
  return MultiFieldFESpace(spaces; style)
end

function evaluate(U::TransientMultiFieldTrialFESpace, t::Real)
  Ut = allocate_trial_space(U)
  evaluate!(Ut, U, t)
  return Ut
end

function evaluate!(Ut::T, U::TransientMultiFieldTrialFESpace, t::Real) where {T}
  spaces_at_t = [evaluate!(Uti, Ui, t) for (Uti, Ui) in zip(Ut, U)]
  style = MultiFieldStyle(U)
  return MultiFieldFESpace(spaces_at_t; style)
end

function evaluate(U::TransientMultiFieldTrialFESpace, t::Nothing)
  spaces = [evaluate(fesp, nothing) for fesp in U.spaces]
  style = MultiFieldStyle(U)
  MultiFieldFESpace(spaces; style)
end

(U::TransientMultiFieldTrialFESpace)(t) = evaluate(U, t)

# Time derivatives
function ∂t(U::TransientMultiFieldTrialFESpace)
  spaces = ∂t.(U.spaces)
  style = MultiFieldStyle(U)
  TransientMultiFieldFESpace(spaces; style)
end

function ∂tt(U::TransientMultiFieldTrialFESpace)
  spaces = ∂tt.(U.spaces)
  style = MultiFieldStyle(U)
  TransientMultiFieldFESpace(spaces; style)
end

# FESpace interface
function zero_free_values(f::TransientMultiFieldTrialFESpace{<:BlockMultiFieldStyle{NB,SB,P}}) where {NB,SB,P}
  block_ranges = get_block_ranges(NB, SB, P)
  block_num_dofs = map(range -> sum(map(num_free_dofs, f.spaces[range])), block_ranges)
  block_vtypes = map(range -> get_vector_type(first(f.spaces[range])), block_ranges)
  return mortar(map(allocate_vector, block_vtypes, block_num_dofs))
end

get_dof_value_type(f::TransientMultiFieldTrialFESpace{MS,CS,V}) where {MS,CS,V} = eltype(V)

get_vector_type(f::TransientMultiFieldTrialFESpace) = f.vector_type

# MultiField interface
# QUESTION Do we need to annotate V?
ConstraintStyle(::Type{TransientMultiFieldTrialFESpace{S,B,V}}) where {S,B,V} = B()
ConstraintStyle(::TransientMultiFieldTrialFESpace) = ConstraintStyle(typeof(f))

# QUESTION Do we need to annotate B and V?
MultiFieldStyle(::Type{TransientMultiFieldTrialFESpace{S,B,V}}) where {S,B,V} = S()
MultiFieldStyle(f::TransientMultiFieldTrialFESpace) = MultiFieldStyle(typeof(f))

function SparseMatrixAssembler(
  mat, vec,
  trial::TransientMultiFieldTrialFESpace{MS},
  test::TransientMultiFieldTrialFESpace{MS},
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

########
# Test #
########
function test_transient_trial_fe_space(Uh)
  UhX = evaluate(Uh, nothing)
  @test UhX isa FESpace

  t = 0.0

  Uh0 = allocate_trial_space(Uh)
  Uh0 = evaluate!(Uh0, Uh, t)
  @test Uh0 isa FESpace

  Uh0 = evaluate(Uh, t)
  @test Uh0 isa FESpace

  Uh0 = Uh(t)
  @test Uh0 isa FESpace

  Uht = ∂t(Uh)
  Uht0 = Uht(t)
  @test Uht0 isa FESpace

  true
end
