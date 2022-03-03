function solve_step!(uf::AbstractVector,
                     solver::ThetaMethod,
                     op::AffineODEOperator,
                     u0::AbstractVector,
                     t0::Real,
                     cache) # -> (uF,tF)

  dt = solver.dt
  solver.θ == 0.0 ? dtθ = dt : dtθ = dt*solver.θ
  tθ = t0+dtθ

  if cache === nothing
    ode_cache = allocate_cache(op)
    vθ = similar(u0)
    vθ .= 0.0
    l_cache = nothing
    A, b = _allocate_matrix_and_vector(op,u0,ode_cache)
  else
    ode_cache, vθ, A, b, l_cache = cache
  end

  ode_cache = update_cache!(ode_cache,op,tθ)

  _matrix_and_vector!(A,b,op,tθ,dtθ,u0,ode_cache,vθ)
  afop = AffineOperator(A,b)

  newmatrix = true
  l_cache = solve!(uf,solver.nls,afop,l_cache,newmatrix)

  uf = uf + u0
  if 0.0 < solver.θ < 1.0
    uf = uf*(1.0/solver.θ)-u0*((1-solver.θ)/solver.θ)
  end

  cache = (ode_cache, vθ, A, b, l_cache)

  tf = t0+dt
  return (uf,tf,cache)

end

function solve_step!(uf::AbstractVector,
                     solver::ThetaMethod,
                     op::ConstantODEOperator,
                     u0::AbstractVector,
                     t0::Real,
                     cache) # -> (uF,tF)

  dt = solver.dt
  solver.θ == 0.0 ? dtθ = dt : dtθ = dt*solver.θ
  tθ = t0+dtθ

  if cache === nothing
    ode_cache = allocate_cache(op)
    vθ = similar(u0)
    vθ .= 0.0
    A, b = _allocate_matrix_and_vector(op,u0,ode_cache)
    A = _matrix!(A,op,tθ,dtθ,u0,ode_cache,vθ)
    b = _vector!(b,op,tθ,dtθ,vθ,ode_cache,vθ)
    M = _allocate_matrix(op,u0,ode_cache)
    M = _mass_matrix!(M,op,tθ,dtθ,u0,ode_cache,vθ)
    _u0 =  similar(u0,(axes(M)[2],)) # Needed for the distributed case
    copy!(_u0,u0)
    l_cache = nothing
    newmatrix = true
  else
    ode_cache, _u0, vθ, A, b, M, l_cache = cache
    newmatrix = false
    copy!(_u0,u0)
  end

  ode_cache = update_cache!(ode_cache,op,tθ)

  vθ = b + M*_u0
  afop = AffineOperator(A,vθ)

  l_cache = solve!(uf,solver.nls,afop,l_cache,newmatrix)

  if 0.0 < solver.θ < 1.0
    uf = uf*(1.0/solver.θ)-u0*((1-solver.θ)/solver.θ)
  end

  cache = (ode_cache, _u0, vθ, A, b, M, l_cache)

  tf = t0+dt
  return (uf,tf,cache)

end

"""
Affine operator that represents the θ-method affine operator at a
given time step, i.e., M(t)(u_n+θ-u_n)/dt + K(t)u_n+θ + b(t)
"""
function ThetaMethodAffineOperator(odeop::AffineODEOperator,tθ::Float64,dtθ::Float64,
                                   u0::AbstractVector,ode_cache,vθ::AbstractVector)
  # vθ .= 0.0
  A, b = _allocate_matrix_and_vector(odeop,u0,ode_cache)
  _matrix_and_vector!(A,b,odeop,tθ,dtθ,u0,ode_cache,vθ)
  afop = AffineOperator(A,b)
end

function _matrix_and_vector!(A,b,odeop,tθ,dtθ,u0,ode_cache,vθ)
  _matrix!(A,odeop,tθ,dtθ,u0,ode_cache,vθ)
  _vector!(b,odeop,tθ,dtθ,u0,ode_cache,vθ)
end

function _matrix!(A,odeop,tθ,dtθ,u0,ode_cache,vθ)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,odeop,tθ,(vθ,vθ),(1.0,1/dtθ),ode_cache)
end

function _mass_matrix!(A,odeop,tθ,dtθ,u0,ode_cache,vθ)
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobian!(A,odeop,tθ,(vθ,vθ),2,(1/dtθ),ode_cache)
end

function _vector!(b,odeop,tθ,dtθ,u0,ode_cache,vθ)
  residual!(b,odeop,tθ,(u0,vθ),ode_cache)
  b .*= -1.0
end

function _allocate_matrix(odeop,u0,ode_cache)
  A = allocate_jacobian(odeop,u0,ode_cache)
  return A
end

function _allocate_matrix_and_vector(odeop,u0,ode_cache)
  b = allocate_residual(odeop,u0,ode_cache)
  A = allocate_jacobian(odeop,u0,ode_cache)
  return A, b
end

"""
Affine operator that represents the θ-method affine operator at a
given time step, i.e., M(t)(u_n+θ-u_n)/dt + K(t)u_n+θ + b(t)
"""
function ThetaMethodConstantOperator(odeop::ConstantODEOperator,tθ::Float64,dtθ::Float64,
                                   u0::AbstractVector,ode_cache,vθ::AbstractVector)
  b = allocate_residual(odeop,u0,ode_cache)
  A = allocate_jacobian(odeop,u0,ode_cache)
  residual!(b,odeop,tθ,(u0,vθ),ode_cache)
  b = -1*b
  z = zero(eltype(A))
  fillstored!(A,z)
  jacobians!(A,odeop,tθ,(vθ,vθ),(1.0,1/dtθ),ode_cache)
  return A, b
end
