function solve_step!(
  x1::NTuple{3,AbstractVector},
  solver::Newmark,
  op::ConstantODEOperator,
  x0::NTuple{3,AbstractVector},
  t0::Real,
  cache) # -> (uF,tF)

  dt = solver.dt
  γ = solver.γ
  β = solver.β
  t1 = t0+dt
  u0, v0, a0 = x0
  u1, v1, a1 = x1

  if cache === nothing
    # Auxiliar variables
    newmatrix = true

    # Allocate caches
    newmark_cache = allocate_cache(op,v0,a0)
    (v,a, ode_cache) = newmark_cache

    # Allocate matrices and vectors
    A, b = _allocate_matrix_and_vector(op,x0,ode_cache)
    M = _allocate_matrix(op,x0,ode_cache)
    C = _allocate_matrix(op,x0,ode_cache)
    b1 = similar(b)
    b1 .= 0.0

    # Define Newmark operator
    newmark_affOp = NewmarkConstantOperator(op,t1,dt,γ,β,(u0,v0,a0),newmark_cache)

    # Fill matrices and vector
    _matrix!(A,newmark_affOp,u1)
    _mass_matrix!(M,newmark_affOp,u1)
    _damping_matrix!(C,newmark_affOp,u1)

    # Create affine operator cache
    affOp_cache = (A,b,M,C,newmark_affOp,nothing)
  else
    newmark_cache, affOp_cache = cache
    newmatrix = false
  end

  # Unpack and update caches
  (v,a, ode_cache) = newmark_cache
  ode_cache = update_cache!(ode_cache,op,t1)
  A,b,M,C,newmark_affOp,l_cache = affOp_cache

  # Update RHS
  _vector!(b,newmark_affOp,u1)
  b1 = b + ( M*(1.0/(β*dt^2)) + C*(γ/(β*dt)) )*u0 +
           ( M*(1.0/(β*dt)) - C*(1-γ/β) )*v0 +
           ( M*(1-2*β)/(2*β) - C*(dt*(1-γ/(2*β))) )*a0

  # Create affine operator with updated RHS
  affOp = AffineOperator(A,b1)
  l_cache = solve!(u1,solver.nls,affOp,l_cache,newmatrix)

  # Update auxiliar variables
  v1 = γ/(β*dt)*(u1-u0) + (1-γ/β)*v0 + dt*(1-γ/(2*β))*a0
  a1 = 1.0/(β*dt^2)*(u1-u0) - 1.0/(β*dt)*v0 - (1-2*β)/(2*β)*a0

  # Pack caches
  affOp_cache = A,b,M,C,newmark_affOp,l_cache
  cache = (newmark_cache, affOp_cache)
  x1 = (u1,v1,a1)

  return (x1,t1,cache)

end

"""
Constant operator that represents the Newmark Affine operator at a
given time step, i.e., M(t)(u_n+1-u_n)/dt + K(t)u_n+1 + b(t)
"""
struct NewmarkConstantOperator <: NonlinearOperator
  odeop::ConstantODEOperator
  t1::Float64
  dt::Float64
  γ::Float64
  β::Float64
  x0::NTuple{3,AbstractVector}
  ode_cache
end

function _matrix_and_vector!(
  A::AbstractMatrix,
  b::AbstractVector,
  affOp::NewmarkConstantOperator,
  x::AbstractVector)
  jacobian!(A,affOp,x)
  residual!(b,affOp,x)
end

function _matrix!(
  A::AbstractMatrix,
  affOp::NewmarkConstantOperator,
  x::AbstractVector)
  jacobian!(A,affOp,x)
end

function _vector!(
  b::AbstractVector,
  affOp::NewmarkConstantOperator,
  x::AbstractVector)
  residual!(b,affOp,x)
end

function residual!(b::AbstractVector,op::NewmarkConstantOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  residual!(b,op.odeop,op.t1,(u1,v1,a1),cache)
  b .*= -1.0
end

function jacobian!(A::AbstractMatrix,op::NewmarkConstantOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.t1,(u1,v1,a1),(1.0,op.γ/(op.β*op.dt),1.0/(op.β*op.dt^2)),cache)
end

function _mass_matrix!(A::AbstractMatrix,op::NewmarkConstantOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobian!(A,op.odeop,op.t1,(u1,v1,a1),3,1.0,cache)
end

function _damping_matrix!(A::AbstractMatrix,op::NewmarkConstantOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobian!(A,op.odeop,op.t1,(u1,v1,a1),2,1.0,cache)
end
