"""
Trait for `ODEOperator` that tells us whether the operator depends on the solution
(including its time derivatives), it is an affine operator that depends on time
or it is a constant operator (affine and time-indepedendent)
"""
abstract type OperatorType end
struct Nonlinear <: OperatorType end
struct Affine  <: OperatorType end
struct Constant  <: OperatorType end
struct ConstantMatrix  <: OperatorType end


"""
It represents the operator in an implicit N-th order ODE, i.e., A(t,u,∂tu,∂t^2u,...,∂t^Nu)
where the implicit PDE reads A(t,u,∂tu,∂t^2u,...,∂t^Nu) = 0, when ∂t^iu is the
i-th time derivative of u, with i=0,..,N. The trait `{C}` determines whether the
operator is fully nonlinear, affine or constant in time.
"""
abstract type ODEOperator{C<:OperatorType} <: GridapType end

"""
It represents an _affine_ operator in an implicit ODE, i.e., an ODE operator of
the form A(t,u,∂tu,...,∂t^Nu) = A_N(t)∂t^Nu + ...A_1(t)∂tu + A_0(t)u + f(t)
"""
const AffineODEOperator = ODEOperator{Affine}

"""
It represents a constant operator in an implicit ODE, i.e., an ODE operator of
the form A(t,u,∂tu,...,∂t^Nu) = A_N∂t^Nu + ...A_1∂tu + A_0u + f
"""
const ConstantODEOperator = ODEOperator{Constant}

"""
It represents an affine operator in an implicit ODE with constant matrix, but
time-dependent right-hand side, i.e., an ODE operator of
the form A(t,u,∂tu,...,∂t^Nu) = A_N∂t^Nu + ...A_1∂tu + A_0u + f(t)
"""
const ConstantMatrixODEOperator = ODEOperator{ConstantMatrix}

"""
Returns the `OperatorType`, i.e., nonlinear, affine, or constant in time
"""
OperatorType(::ODEOperator{C}) where {C} = C

"""
Returns the order of the ODE operator
"""
function get_order(::ODEOperator)
  @abstractmethod
end

"""
It provides A(t,u,∂tu,...,∂t^Nu) for a given (t,u,∂tu,...,∂t^Nu)
"""
function residual!(
  r::AbstractVector,
  op::ODEOperator,
  t::Real,
  u::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  ode_cache)
  @abstractmethod
end

"""
"""
function allocate_residual(
  op::ODEOperator,
  u::Union{AbstractVector,Tuple{Vararg{AbstractVector}}},
  ode_cache)
  @abstractmethod
end

"""
It adds contribution to the Jacobian with respect to the i-th time derivative,
with i=0,...,N. That is, adding γ_i*[∂A/∂(∂t^iu)](t,u,∂tu,...,∂t^Nu) for a
given (t,u,∂tu,...,∂t^Nu) to a given matrix J, where γ_i is a scaling coefficient
provided by the `ODESolver`, e.g., 1/Δt for Backward Euler; It represents
∂(δt^i(u))/∂(u), in which δt^i(⋅) is the approximation of ∂t^i(⋅) in the solver.
Note that for i=0, γ_i=1.0.
"""
function jacobian!(
  J::AbstractMatrix,
  op::ODEOperator,
  t::Real,
  u::Tuple{Vararg{AbstractVector}},
  i::Integer,
  γᵢ::Real,
  ode_cache)
  @abstractmethod
  # Add values to J
end

"""
Add the contribution of all jacobians ,i.e., ∑ᵢ γ_i*[∂A/∂(∂t^iu)](t,u,∂tu,...,∂t^Nu)
"""
function jacobians!(
  J::AbstractMatrix,
  op::ODEOperator,
  t::Real,
  u::Tuple{Vararg{AbstractVector}},
  γ::Tuple{Vararg{Real}},
  ode_cache)
  @abstractmethod
  # Add values to J
end

"""
"""
function allocate_jacobian(op::ODEOperator,u::AbstractVector,ode_cache)
  @abstractmethod
end

"""
Allocates the cache data required by the `ODESolution` for a given `ODEOperator`
"""
allocate_cache(op::ODEOperator) = @abstractmethod

#@fverdugo to be used as `cache = update_cache!(cache,op,t)`
update_cache!(cache,op::ODEOperator,t::Real) = @abstractmethod

"""
Tests the interface of `ODEOperator` specializations
"""
function test_ode_operator(op::ODEOperator,t::Real,u::AbstractVector,u_t::AbstractVector)
  cache = allocate_cache(op)
  cache = update_cache!(cache,op,0.0)
  r = allocate_residual(op,u,cache)
  residual!(r,op,t,(u,u_t),cache)
  J = allocate_jacobian(op,u,cache)
  jacobian!(J,op,t,(u,u_t),1,1.0,cache)
  jacobian!(J,op,t,(u,u_t),2,1.0,cache)
  jacobians!(J,op,t,(u,u_t),(1.0,1.0),cache)
  true
end
