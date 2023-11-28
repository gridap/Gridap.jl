using SparseArrays: spzeros

using Gridap
using Gridap.Algebra
using Gridap.Polynomials
using Gridap.ODEs

# Toy 1st-order linear ODE
# u_1_t - a * u_1 = 0
# u_2_t - b * u_1 - c * u_2 = 0
# The analytical solution is
# If a = c
# u_1(t) = A exp(a * t)
# u_2(t) = (B + b A t) exp(a * t)
# If a != c
# u_1(t) = A exp(a * t)
# u_2(t) = b / (a - c) A exp(a * t) + B * exp(c * t)

# Toy 2nd-order linear ODE
# u_1_tt + b * u_1_t - a * u_1 = 0
# u_2_tt + a * u_1_t - b * u_1 - c * u_2 = 0

struct ODEOperatorMock{T,C} <: ODEOperator{C}
  a::T
  b::T
  c::T
  order::Integer

  function ODEOperatorMock{C}(a, b, c, order) where {C}
    T = promote_type(typeof(a), typeof(b), typeof(c))
    new{T,C}(a, b, c, order)
  end
end

Gridap.Polynomials.get_order(op::ODEOperatorMock) = op.order

function Gridap.Algebra.allocate_residual(
  op::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  zeros(2)
end

function Gridap.Algebra.residual!(
  r::AbstractVector, op::ODEOperatorMock,
  t::Real, us::NTuple{2,AbstractVector},
  cache
)
  u, u̇ = us
  a, b, c = op.a, op.b, op.c
  r[1] = u̇[1] - a * u[1]
  r[2] = u̇[2] - b * u[1] - c * u[2]
  r
end

function Gridap.Algebra.residual!(
  r::AbstractVector, op::ODEOperatorMock,
  t::Real, us::NTuple{3,AbstractVector},
  cache
)
  u, u̇, ü = us
  a, b, c = op.a, op.b, op.c
  r[1] = ü[1] + b * u̇[1] - a * u[1]
  r[2] = ü[2] + a * u̇[1] - b * u[1] - c * u[2]
  r
end

function Gridap.Algebra.allocate_jacobian(
  op::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  spzeros(2, 2)
end

function Gridap.Algebra.jacobian!(
  J::AbstractMatrix, op::ODEOperatorMock,
  t::Real, us::NTuple{2,AbstractVector},
  k::Integer, γ::Real,
  cache
)
  @assert get_order(op) == 1
  @assert 0 <= k < get_order(op) + 1

  a, b, c = op.a, op.b, op.c
  if k == 0
    J[1, 1] += γ * (-a)
    J[2, 1] += γ * (-b)
    J[2, 2] += γ * (-c)
  elseif k == 1
    J[1, 1] += γ
    J[2, 2] += γ
  end
  J
end

function Gridap.Algebra.jacobian!(
  J::AbstractMatrix, op::ODEOperatorMock,
  t::Real, us::NTuple{3,AbstractVector},
  k::Integer, γ::Real,
  cache
)
  @assert get_order(op) == 2
  @assert 0 <= k < get_order(op) + 1

  a, b, c = op.a, op.b, op.c
  if k == 0
    J[1, 1] += γ * (-a)
    J[2, 1] += γ * (-b)
    J[2, 2] += γ * (-c)
  elseif k == 1
    J[1, 1] += γ * b
    J[2, 2] += γ * a
  elseif k == 2
    J[1, 1] += γ
    J[2, 2] += γ
  end
  J
end
