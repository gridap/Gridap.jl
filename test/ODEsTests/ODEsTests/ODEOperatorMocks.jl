# Toy linear ODE with 2 DOFs
# u_1_t - a * u_1 = 0
# u_2_t - b * u_1 - c * u_2 = 0

# Toy 2nd order ODE with 2 DOFs
# u_1_tt + b * u_1_t - a * u_1 = 0
# u_2_tt + a * u_1_t - b * u_1 - c * u_2 = 0

import Gridap.ODEs.ODETools: ODEOperator
import Gridap.ODEs.ODETools: AffineODEOperator
import Gridap.ODEs.ODETools: ConstantODEOperator
import Gridap.ODEs.ODETools: allocate_cache
import Gridap.ODEs.ODETools: update_cache!
import Gridap.ODEs.ODETools: allocate_residual
import Gridap.ODEs.ODETools: jacobian!
import Gridap.ODEs.ODETools: jacobians!
import Gridap.ODEs.ODETools: allocate_jacobian
import Gridap.ODEs.ODETools: residual!
import Gridap.ODEs.ODETools: rhs!
import Gridap.ODEs.ODETools: explicit_rhs!
import Gridap.ODEs.ODETools: lhs!
using SparseArrays: spzeros

struct ODEOperatorMock{T<:Real,C} <: ODEOperator{C}
  a::T
  b::T
  c::T
  order::Integer
end

get_order(op::ODEOperatorMock) = op.order

function residual!(r::AbstractVector,op::ODEOperatorMock,t::Real,x::NTuple{2,AbstractVector},ode_cache)
  u,u_t = x
  r .= 0
  r[1] = u_t[1] - op.a * u[1]
  r[2] = u_t[2] - op.b * u[1] - op.c * u[2]
  r
end

function rhs!(r::AbstractVector,op::ODEOperatorMock,t::Real,x::NTuple{2,AbstractVector},ode_cache)
  u,u_t = x
  r .= 0
  r[1] = op.a * u[1]
  r[2] = op.b * u[1] + op.c * u[2]
  r
end

function explicit_rhs!(r::AbstractVector,op::ODEOperatorMock,t::Real,x::NTuple{2,AbstractVector},ode_cache)
  u,u_t = x
  r .= 0
  r
end

function lhs!(r::AbstractVector,op::ODEOperatorMock,t::Real,x::NTuple{2,AbstractVector},ode_cache)
  u,u_t = x
  r .= 0
  r[1] = u_t[1]
  r[2] = u_t[2]
  r
end

function residual!(r::AbstractVector,op::ODEOperatorMock,t::Real,x::NTuple{3,AbstractVector},ode_cache)
  u,u_t,u_tt = x
  r .= 0
  r[1] = u_tt[1] + op.b * u_t[1] - op.a * u[1]
  r[2] = u_tt[2] + op.a * u_t[1]- op.b * u[1] - op.c * u[2]
  r
end

function allocate_residual(op::ODEOperatorMock,t0::Real,u::AbstractVector,cache)
  zeros(2)
end

function jacobian!(J::AbstractMatrix,
  op::ODEOperatorMock,
  t::Real,
  x::NTuple{2,AbstractVector},
  i::Int,
  γᵢ::Real,
  ode_cache)
  @assert get_order(op) == 1
  @assert 0 < i <= get_order(op)+1
  if i==1
    J[1,1] += -op.a*γᵢ
    J[2,1] += -op.b*γᵢ
    J[2,2] += -op.c*γᵢ
  elseif i==2
    J[1,1] += 1.0*γᵢ
    J[2,2] += 1.0*γᵢ
  end
  J
end

function jacobian!(J::AbstractMatrix,
  op::ODEOperatorMock,
  t::Real,
  x::NTuple{3,AbstractVector},
  i::Int,
  γᵢ::Real,
  ode_cache)
  @assert get_order(op) == 2
  @assert 0 < i <= get_order(op)+1
  if i==1
    J[1,1] += -op.a*γᵢ
    J[2,1] += -op.b*γᵢ
    J[2,2] += -op.c*γᵢ
  elseif i==2
    J[1,1] += op.b*γᵢ
    J[2,2] += op.a*γᵢ
  elseif i==3
    J[1,1] += 1.0*γᵢ
    J[2,2] += 1.0*γᵢ
  end
  J
end

function jacobians!(
  J::AbstractMatrix,
  op::ODEOperatorMock,
  t::Real,
  x::Tuple{Vararg{AbstractVector}},
  γ::Tuple{Vararg{Real}},
  ode_cache)
  @assert length(γ) == get_order(op) + 1
  for order in 1:get_order(op)+1
    jacobian!(J,op,t,x,order,γ[order],ode_cache)
  end
  J
end

function allocate_jacobian(op::ODEOperatorMock,t0::Real,u::AbstractVector,cache)
  spzeros(2,2)
end

allocate_cache(op::ODEOperatorMock) = nothing
allocate_cache(op::ODEOperatorMock,v::AbstractVector) = (similar(v),nothing)
allocate_cache(op::ODEOperatorMock,v::AbstractVector,a::AbstractVector) = (similar(v),similar(a),nothing)
update_cache!(cache,op::ODEOperatorMock,t::Real) = cache
