using SparseArrays: spzeros

using Gridap
using Gridap.Algebra
using Gridap.Polynomials
using Gridap.ODEs

###################
# ODEOperatorMock #
###################
"""
    struct ODEOperatorMock <: ODEOperator end

Mock linear ODE of arbitrary order
```math
∑_{0 ≤ k ≤ N} form_k(t) ∂t^k u + forcing(t) = 0.
```
"""
struct ODEOperatorMock{T} <: ODEOperator{T}
  forms::Tuple{Vararg{Function}}
  forcing::Function
end

Polynomials.get_order(odeop::ODEOperatorMock) = length(odeop.forms) - 1

function Algebra.allocate_residual(
  odeop::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  f = forcing(t)
  copy(f)
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  fill!(r, zero(eltype(r)))
  r .+= odeop.forcing(t)
  for k in 0:get_order(odeop)
    mat = odeop.forms[k+1](t)
    r .+= mat * us[k+1]
  end
  r
end

function Algebra.allocate_jacobian(
  odeop::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  T = eltype(first(us))
  n = length(first(us))
  spzeros(T, n, n)
end

function Algebra.jacobian!(
  J::AbstractMatrix, odeop::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  k::Integer, γ::Real,
  odeopcache
)
  @assert 0 <= k <= get_order(odeop)
  jac = odeop.forms[k+1](t)
  @. J += γ * jac
  J
end
