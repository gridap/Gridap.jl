using SparseArrays: spzeros

using Gridap
using Gridap.Algebra
using Gridap.Polynomials
using Gridap.ODEs

# Toy 1st-order linear ODE
# M u̇ + K u + f(t) = 0
# u(t) = exp(-M\K t) [u0 - ∫_{0, t} exp(M\K s) M\f(s) ds

struct ODEOperatorMock{C} <: ODEOperator{C}
  M::AbstractMatrix
  K::AbstractMatrix
  f::Function
end

Polynomials.get_order(op::ODEOperatorMock) = 1

function Algebra.allocate_residual(
  op::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  M = op.M
  zeros(typeof(t), size(M, 2))
end

function Algebra.residual!(
  r::AbstractVector, op::ODEOperatorMock,
  t::Real, us::NTuple{2,AbstractVector},
  cache; include_highest::Bool=true
)
  u, u̇ = us
  fill!(r, zero(eltype(r)))
  if include_highest
    r .+= op.M * u̇
  end
  r .+= op.K * u
  r .+= f(t)
  r
end

function Algebra.allocate_jacobian(
  op::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  M = op.M
  spzeros(typeof(t), size(M, 2), size(M, 2))
end

function Algebra.jacobian!(
  J::AbstractMatrix, op::ODEOperatorMock,
  t::Real, us::NTuple{2,AbstractVector},
  k::Integer, γ::Real,
  cache
)
  if k == 0
    @. J += γ * op.K
  elseif k == 1
    @. J += γ * op.M
  end
  J
end
