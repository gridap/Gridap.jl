function solve_step!(
  x1::NTuple{3,AbstractVector},
  solver::Newmark,
  op::ConstantMatrixODEOperator,
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
    newmatrix = true

    # Allocate caches
    newmark_cache = allocate_cache(op,v0,a0)
    (v,a, ode_cache) = newmark_cache

    # Define Newmark operator
    newmark_affOp = NewmarkConstantMatrixOperator(op,t1,dt,γ,β,(u0,v0,a0),newmark_cache)

    # Allocate matrices and vectors
    A, b = _allocate_matrix_and_vector(op,x0,ode_cache)
    jacobian!(A,newmark_affOp,u1)

    # Create affine operator cache
    affOp_cache = (A,b,newmark_affOp,nothing)
  else
    newmatrix = false
    newmark_cache, affOp_cache = cache
  end

  # Unpack and update caches
  (v,a, ode_cache) = newmark_cache
  ode_cache = update_cache!(ode_cache,op,t1)
  A,b,newmark_affOp,l_cache = affOp_cache

  # Fill vector
  newmark_affOp = NewmarkConstantMatrixOperator(op,t1,dt,γ,β,x0,newmark_cache)
  residual!(b,newmark_affOp,u1)

  # Create affine operator with updated RHS
  affOp = AffineOperator(A,b)
  l_cache = solve!(u1,solver.nls,affOp,l_cache,newmatrix)

  # Update auxiliar variables
  u1 = u1 + u0
  v1 = γ/(β*dt)*(u1-u0) + (1-γ/β)*v0 + dt*(1-γ/(2*β))*a0
  a1 = 1.0/(β*dt^2)*(u1-u0) - 1.0/(β*dt)*v0 - (1-2*β)/(2*β)*a0

  # Pack caches
  affOp_cache = A,b,newmark_affOp,l_cache
  cache = (newmark_cache, affOp_cache)
  x1 = (u1,v1,a1)

  return (x1,t1,cache)

end

"""
Affine operator that represents the Newmark Affine operator with constant
matrix at a given time step, i.e., M(u_n+1-u_n)/dt + K u_n+1 + b(t)
"""
mutable struct NewmarkConstantMatrixOperator <: NonlinearOperator
  odeop::ConstantMatrixODEOperator
  t1::Float64
  dt::Float64
  γ::Float64
  β::Float64
  x0::NTuple{3,AbstractVector}
  ode_cache
end

function residual!(b::AbstractVector,op::NewmarkConstantMatrixOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  residual!(b,op.odeop,op.t1,(u1,v1,a1),cache)
  b .*= -1.0
end

function jacobian!(A::AbstractMatrix,op::NewmarkConstantMatrixOperator,x::AbstractVector)
  u1 = x
  u0, v0, a0 = op.x0
  v1, a1, cache = op.ode_cache
  a1 = 1.0/(op.β*op.dt^2)*(u1-u0) - 1.0/(op.β*op.dt)*v0 - (1-2*op.β)/(2*op.β)*a0
  v1 = op.γ/(op.β*op.dt)*(u1-u0) + (1-op.γ/op.β)*v0 + op.dt*(1-op.γ/(2*op.β))*a0
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,op.odeop,op.t1,(u1,v1,a1),(1.0,op.γ/(op.β*op.dt),1.0/(op.β*op.dt^2)),cache)
end

function _allocate_matrix(odeop::NewmarkConstantMatrixOperator,x::Tuple{Vararg{AbstractVector}},ode_cache)
  A = allocate_jacobian(odeop,x[1],ode_cache)
  return A
end

function _allocate_matrix_and_vector(odeop::NewmarkConstantMatrixOperator,x::Tuple{Vararg{AbstractVector}},ode_cache)
  b = allocate_residual(odeop,x[1],ode_cache)
  A = allocate_jacobian(odeop,x[1],ode_cache)
  return A, b
end
