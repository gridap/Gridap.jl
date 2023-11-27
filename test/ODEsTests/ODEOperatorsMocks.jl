using SparseArrays: spzeros

using Gridap
using Gridap.Algebra
using Gridap.ODEs

import Gridap.Algebra: allocate_residual
import Gridap.Algebra: residual!
import Gridap.Algebra: allocate_jacobian
import Gridap.Algebra: jacobian!

import Gridap.ODEs: ODEOperator
import Gridap.ODEs: allocate_cache
import Gridap.ODEs: update_cache!
import Gridap.ODEs: jacobians!

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
  order::Int

  function ODEOperatorMock{C}(a, b, c, order) where {C}
    T = promote_type(typeof(a), typeof(b), typeof(c))
    new{T,C}(a, b, c, order)
  end
end

get_order(op::ODEOperatorMock) = op.order

function allocate_cache(op::ODEOperatorMock)
  nothing
end

function update_cache!(cache, op::ODEOperatorMock, t::Real)
  cache
end

function allocate_residual(
  op::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  zeros(2)
end

function residual!(
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

function residual!(
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

function allocate_jacobian(
  op::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  spzeros(2, 2)
end

function jacobian!(
  J::AbstractMatrix, op::ODEOperatorMock,
  t::Real, us::NTuple{2,AbstractVector},
  i::Int, γ::Real,
  cache
)
  @assert get_order(op) == 1
  @assert 0 <= i < get_order(op) + 1

  a, b, c = op.a, op.b, op.c
  if i == 0
    J[1, 1] += γ * (-a)
    J[2, 1] += γ * (-b)
    J[2, 2] += γ * (-c)
  elseif i == 1
    J[1, 1] += γ
    J[2, 2] += γ
  end
  J
end

function jacobian!(
  J::AbstractMatrix, op::ODEOperatorMock,
  t::Real, us::NTuple{3,AbstractVector},
  i::Int, γ::Real,
  cache
)
  @assert get_order(op) == 2
  @assert 0 <= i < get_order(op) + 1

  a, b, c = op.a, op.b, op.c
  if i == 0
    J[1, 1] += γ * (-a)
    J[2, 1] += γ * (-b)
    J[2, 2] += γ * (-c)
  elseif i == 1
    J[1, 1] += γ * b
    J[2, 2] += γ * a
  elseif i == 2
    J[1, 1] += γ
    J[2, 2] += γ
  end
  J
end

function jacobians!(
  J::AbstractMatrix, op::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  cache
)
  @assert length(γs) == get_order(op) + 1
  for order in 0:get_order(op)
    γ = γs[order+1]
    jacobian!(J, op, t, us, order, γ, cache)
  end
  J
end
