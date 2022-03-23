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
  Ut
end

"""
We can evaluate at `nothing` when we do not care about the Dirichlet vals
"""
function evaluate(U::TransientTrialFESpace,t::Nothing)
  U.Ud0
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

function evaluate!(Ut::FESpace,U::FESpace,t::Real)
  U
end

function allocate_trial_space(U::FESpace)
  U
end

function evaluate(U::FESpace,t::Real)
  U
end

function evaluate(U::FESpace,t::Nothing)
  U
end

@static if VERSION >= v"1.3"
  (U::FESpace)(t) = U
end

# Define the interface for MultiField

struct TransientMultiFieldTrialFESpace
  spaces::Vector
end
Base.iterate(m::TransientMultiFieldTrialFESpace) = iterate(m.spaces)
Base.iterate(m::TransientMultiFieldTrialFESpace,state) = iterate(m.spaces,state)
Base.getindex(m::TransientMultiFieldTrialFESpace,field_id::Integer) = m.spaces[field_id]
Base.length(m::TransientMultiFieldTrialFESpace) = length(m.spaces)


function TransientMultiFieldFESpace(spaces::Vector)
  TransientMultiFieldTrialFESpace(spaces)
end

function TransientMultiFieldFESpace(spaces::Vector{<:SingleFieldFESpace})
  MultiFieldFESpace(spaces)
end

function evaluate!(Ut::T,U::TransientMultiFieldTrialFESpace,t::Real) where T
  spaces_at_t = [evaluate!(Uti,Ui,t) for (Uti,Ui) in zip(Ut,U)]
  MultiFieldFESpace(spaces_at_t)
end

function allocate_trial_space(U::TransientMultiFieldTrialFESpace)
  spaces = allocate_trial_space.(U.spaces)
  MultiFieldFESpace(spaces)
end

function evaluate(U::TransientMultiFieldTrialFESpace,t::Real)
  Ut = allocate_trial_space(U)
  evaluate!(Ut,U,t)
  Ut
end

function evaluate(U::TransientMultiFieldTrialFESpace,t::Nothing)
  MultiFieldFESpace([evaluate(fesp,nothing) for fesp in U.spaces])
end

(U::TransientMultiFieldTrialFESpace)(t) = evaluate(U,t)

function ∂t(U::TransientMultiFieldTrialFESpace)
  spaces = ∂t.(U.spaces)
  TransientMultiFieldFESpace(spaces)
end
