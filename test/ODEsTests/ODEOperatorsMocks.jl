using LinearAlgebra
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
```
∑_{0 ≤ k ≤ N} form_k(t) ∂t^k u - forcing(t) = 0.
```
"""
struct ODEOperatorMock{T} <: ODEOperator{T}
  forms::Tuple{Vararg{Function}}
  forcing::Function
end

Polynomials.get_order(odeop::ODEOperatorMock) = length(odeop.forms) - 1

ODEs.get_forms(odeop::ODEOperatorMock{<:AbstractQuasilinearODE}) = (odeop.forms[end],)
ODEs.get_forms(odeop::ODEOperatorMock{<:AbstractLinearODE}) = odeop.forms

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
  odeopcache; add::Bool=false
)
  order = get_order(odeop)
  !add && fill!(r, zero(eltype(r)))
  axpy!(-1, odeop.forcing(t), r)
  for k in 0:order
    mat = odeop.forms[k+1](t)
    axpy!(1, mat * us[k+1], r)
  end
  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOperatorMock{<:AbstractQuasilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  order = get_order(odeop)
  !add && fill!(r, zero(eltype(r)))
  axpy!(-1, odeop.forcing(t), r)
  for k in 0:order-1
    mat = odeop.forms[k+1](t)
    axpy!(1, mat * us[k+1], r)
  end
  k = order
  mat = odeop.forms[k+1](t)
  axpy!(1, mat * us[k+1], r)
  r
end

function Algebra.residual!(
  r::AbstractVector, odeop::ODEOperatorMock{<:AbstractLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  order = get_order(odeop)
  !add && fill!(r, zero(eltype(r)))
  axpy!(-1, odeop.forcing(t), r)
  for k in 0:order
    mat = odeop.forms[k+1](t)
    axpy!(1, mat * us[k+1], r)
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
  J = spzeros(T, n, n)
  fill!(J, 1)
  J
end

function ODEs.jacobian_add!(
  J::AbstractMatrix, odeop::ODEOperatorMock,
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  for (w, form) in zip(ws, odeop.forms)
    iszero(w) && continue
    jac = form(t)
    axpy_entries!(w, jac, J)
  end
end
