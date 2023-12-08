"""
A single field FE space with transient Dirichlet data (see Multifield below).
"""
struct TransientTrialFESpace{A,B}
  space::A
  dirichlet_t::Union{Function,Vector{<:Function}}
  Ud0::B

  function TransientTrialFESpace(space::A,dirichlet_t::Union{Function,Vector{<:Function}}) where A
    Ud0 = HomogeneousTrialFESpace(space)
    B = typeof(Ud0)
    new{A,B}(space,dirichlet_t,Ud0)
  end
end

function TransientTrialFESpace(space::A) where A
  HomogeneousTrialFESpace(space)
end

"""
Time evaluation without allocating Dirichlet vals
"""
function evaluate!(Ut::T,U::TransientTrialFESpace,t::Real) where T
  if isa(U.dirichlet_t,Vector)
    objects_at_t = map( o->o(t), U.dirichlet_t)
  else
    objects_at_t = U.dirichlet_t(t)
  end
  TrialFESpace!(Ut,objects_at_t)
  Ut
end

"""
Allocate the space to be used as first argument in evaluate!
"""
function allocate_trial_space(U::TransientTrialFESpace)
  HomogeneousTrialFESpace(U.space)
end

"""
Time evaluation allocating Dirichlet vals
"""
function evaluate(U::TransientTrialFESpace,t::Real)
  Ut = allocate_trial_space(U)
  evaluate!(Ut,U,t)
  return Ut
end

"""
We can evaluate at `nothing` when we do not care about the Dirichlet vals
"""
function evaluate(U::TransientTrialFESpace,t::Nothing)
  return U.Ud0
end

evaluate(U::TrialFESpace,t::Nothing) = U

"""
Functor-like evaluation. It allocates Dirichlet vals in general.
"""
(U::TransientTrialFESpace)(t) = evaluate(U,t)

(U::TrialFESpace)(t) = U
(U::ZeroMeanFESpace)(t) = U
# (U::Union{TrialFESpace,ZeroMeanFESpace})(t) = U

"""
Time derivative of the Dirichlet functions
"""
∂t(U::TransientTrialFESpace) = TransientTrialFESpace(U.space,∂t.(U.dirichlet_t))
∂t(U::SingleFieldFESpace) = HomogeneousTrialFESpace(U)
∂t(U::MultiFieldFESpace) = MultiFieldFESpace(∂t.(U.spaces))
∂t(t::T) where T<:Number = zero(T)

"""
Time 2nd derivative of the Dirichlet functions
"""
∂tt(U::TransientTrialFESpace) = TransientTrialFESpace(U.space,∂tt.(U.dirichlet_t))
∂tt(U::SingleFieldFESpace) = HomogeneousTrialFESpace(U)
∂tt(U::MultiFieldFESpace) = MultiFieldFESpace(∂tt.(U.spaces))
∂tt(t::T) where T<:Number = zero(T)

zero_free_values(f::TransientTrialFESpace) = zero_free_values(f.space)
has_constraints(f::TransientTrialFESpace) = has_constraints(f.space)
get_dof_value_type(f::TransientTrialFESpace) = get_dof_value_type(f.space)
get_vector_type(f::TransientTrialFESpace) = get_vector_type(f.space)

# Testing the interface

function test_transient_trial_fe_space(Uh)
  UhX = evaluate(Uh,nothing)
  @test isa(UhX,FESpace)
  Uh0 = allocate_trial_space(Uh)
  Uh0 = evaluate!(Uh0,Uh,0.0)
  @test isa(Uh0,FESpace)
  Uh0 = evaluate(Uh,0.0)
  @test isa(Uh0,FESpace)
  Uh0 = Uh(0.0)
  @test isa(Uh0,FESpace)
  Uht=∂t(Uh)
  Uht0=Uht(0.0)
  @test isa(Uht0,FESpace)
  true
end

# Define the TransientTrialFESpace interface for stationary spaces

evaluate!(Ut::FESpace,U::FESpace,t::Real) = U
allocate_trial_space(U::FESpace) = U
evaluate(U::FESpace,t::Real) = U
evaluate(U::FESpace,t::Nothing) = U

@static if VERSION >= v"1.3"
  (U::FESpace)(t) = U
end

# Define the interface for MultiField

struct TransientMultiFieldTrialFESpace{MS<:MultiFieldStyle,CS<:ConstraintStyle,V}
  vector_type::Type{V}
  spaces::Vector
  multi_field_style::MS
  constraint_style::CS
  function TransientMultiFieldTrialFESpace(
    ::Type{V},
    spaces::Vector,
    multi_field_style::MultiFieldStyle) where V
    @assert length(spaces) > 0

    MS = typeof(multi_field_style)
    if any( map(has_constraints,spaces) )
      constraint_style = Constrained()
    else
      constraint_style = UnConstrained()
    end
    CS = typeof(constraint_style)
    new{MS,CS,V}(V,spaces,multi_field_style,constraint_style)
  end
end

# Default constructors
function TransientMultiFieldFESpace(spaces::Vector; 
                                    style = ConsecutiveMultiFieldStyle())
  Ts = map(get_dof_value_type,spaces)
  T  = typeof(*(map(zero,Ts)...))
  if isa(style,BlockMultiFieldStyle)
    style = BlockMultiFieldStyle(style,spaces)
    VT = typeof(mortar(map(zero_free_values,spaces)))
  else
    VT = Vector{T}
  end
  TransientMultiFieldTrialFESpace(VT,spaces,style)
end

function TransientMultiFieldFESpace(::Type{V},spaces::Vector) where V
  TransientMultiFieldTrialFESpace(V,spaces,ConsecutiveMultiFieldStyle())
end

function TransientMultiFieldFESpace(spaces::Vector{<:SingleFieldFESpace}; 
                                    style = ConsecutiveMultiFieldStyle())
  MultiFieldFESpace(spaces,style=style)
end

function TransientMultiFieldFESpace(::Type{V},spaces::Vector{<:SingleFieldFESpace}) where V
  MultiFieldFESpace(V,spaces,ConsecutiveMultiFieldStyle())
end

Base.iterate(m::TransientMultiFieldTrialFESpace) = iterate(m.spaces)
Base.iterate(m::TransientMultiFieldTrialFESpace,state) = iterate(m.spaces,state)
Base.getindex(m::TransientMultiFieldTrialFESpace,field_id::Integer) = m.spaces[field_id]
Base.length(m::TransientMultiFieldTrialFESpace) = length(m.spaces)

function evaluate!(Ut::T,U::TransientMultiFieldTrialFESpace,t::Real) where T
  spaces_at_t = [evaluate!(Uti,Ui,t) for (Uti,Ui) in zip(Ut,U)]
  mfs = MultiFieldStyle(U)
  return MultiFieldFESpace(spaces_at_t;style=mfs)
end

function allocate_trial_space(U::TransientMultiFieldTrialFESpace)
  spaces = allocate_trial_space.(U.spaces)
  mfs = MultiFieldStyle(U)
  return MultiFieldFESpace(spaces;style=mfs)
end

function evaluate(U::TransientMultiFieldTrialFESpace,t::Real)
  Ut = allocate_trial_space(U)
  evaluate!(Ut,U,t)
  return Ut
end

function evaluate(U::TransientMultiFieldTrialFESpace,t::Nothing)
  spaces = [evaluate(fesp,nothing) for fesp in U.spaces]
  mfs = MultiFieldStyle(U)
  MultiFieldFESpace(spaces;style=mfs)
end

(U::TransientMultiFieldTrialFESpace)(t) = evaluate(U,t)

function ∂t(U::TransientMultiFieldTrialFESpace)
  spaces = ∂t.(U.spaces)
  mfs = MultiFieldStyle(U)
  TransientMultiFieldFESpace(spaces;style=mfs)
end

function zero_free_values(f::TransientMultiFieldTrialFESpace{<:BlockMultiFieldStyle{NB,SB,P}}) where {NB,SB,P}
  block_ranges   = get_block_ranges(NB,SB,P)
  block_num_dofs = map(range->sum(map(num_free_dofs,f.spaces[range])),block_ranges)
  block_vtypes   = map(range->get_vector_type(first(f.spaces[range])),block_ranges)
  values = mortar(map(allocate_vector,block_vtypes,block_num_dofs))
  fill!(values,zero(eltype(values)))
  return values
end

get_dof_value_type(f::TransientMultiFieldTrialFESpace{MS,CS,V}) where {MS,CS,V} = eltype(V)
get_vector_type(f::TransientMultiFieldTrialFESpace) = f.vector_type
ConstraintStyle(::Type{TransientMultiFieldTrialFESpace{S,B,V}}) where {S,B,V} = B()
ConstraintStyle(::TransientMultiFieldTrialFESpace) = ConstraintStyle(typeof(f))
MultiFieldStyle(::Type{TransientMultiFieldTrialFESpace{S,B,V}}) where {S,B,V} = S()
MultiFieldStyle(f::TransientMultiFieldTrialFESpace) = MultiFieldStyle(typeof(f))

function SparseMatrixAssembler(mat,vec,
                               trial::TransientMultiFieldTrialFESpace{MS},
                               test ::TransientMultiFieldTrialFESpace{MS},
                               strategy::AssemblyStrategy=DefaultAssemblyStrategy()
                               ) where MS <: BlockMultiFieldStyle
  mfs = MultiFieldStyle(test)
  return BlockSparseMatrixAssembler(mfs,trial,test,SparseMatrixBuilder(mat),ArrayBuilder(vec),strategy)
end
