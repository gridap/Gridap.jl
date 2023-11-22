"""
A wrapper of `TransientFEOperator` that transforms it to `ODEOperator`, i.e.,
takes A(t,uh,∂tuh,∂t^2uh,...,∂t^Nuh,vh) and returns A(t,uF,∂tuF,...,∂t^NuF)
where uF,∂tuF,...,∂t^NuF represent the free values of the `EvaluationFunction`
uh,∂tuh,∂t^2uh,...,∂t^Nuh.
"""
struct ODEOpFromFEOp{C} <: ODEOperator{C}
  feop::TransientFEOperator{C}
end

get_order(op::ODEOpFromFEOp) = get_order(op.feop)

# ODEOperator interface
function allocate_residual(op::ODEOpFromFEOp, t::Real, u::AbstractVector, cache)
  Us, Uts, fecache = cache
  uh = EvaluationFunction(Us[1], u)
  allocate_residual(op.feop, t, uh, fecache)
end

function residual!(
  r::AbstractVector, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  cache
)
  uh = _make_uh_from_us(op, us, cache)
  residual!(r, op.feop, t, uh, cache)
end

function allocate_jacobian(op::ODEOpFromFEOp, t::Real, u::AbstractVector, cache)
  Us, Uts, fecache = cache
  uh = EvaluationFunction(Us[1], u)
  allocate_jacobian(op.feop, t, uh, fecache)
end

function jacobian!(
  J::AbstractMatrix, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  i::Integer, γ::Real,
  cache
)
  uh = _make_uh_from_us(op, us, cache)
  jacobian!(J, op.feop, t, uh, i, γ, cache)
end

function jacobians!(
  J::AbstractMatrix, op::ODEOpFromFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  γs::Tuple{Vararg{Real}},
  cache
)
  uh = _make_uh_from_us(op, us, cache)
  jacobians!(J, op.feop, t, uh, γs, cache)
end

function allocate_cache(op::ODEOpFromFEOp)
  Ut = get_trial(op.feop)
  U = allocate_trial_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for i in 1:get_order(op)
    Uts = (Uts..., ∂t(Uts[i]))
    Us = (Us..., allocate_trial_space(Uts[i+1]))
  end
  fecache = allocate_cache(op.feop)
  ode_cache = (Us, Uts, fecache)
  ode_cache
end

function allocate_cache(op::ODEOpFromFEOp, v::AbstractVector)
  ode_cache = allocate_cache(op)
  _v = similar(v)
  (_v, ode_cache)
end

function allocate_cache(op::ODEOpFromFEOp, v::AbstractVector, a::AbstractVector)
  ode_cache = allocate_cache(op)
  _v = similar(v)
  _a = similar(a)
  (_v, _a, ode_cache)
end

function update_cache!(ode_cache, op::ODEOpFromFEOp, t::Real)
  _Us, Uts, fecache = ode_cache
  Us = ()
  for i in 1:get_order(op)+1
    Us = (Us..., evaluate!(_Us[i], Uts[i], t))
  end
  fecache = update_cache!(fecache, op.feop, t)
  (Us, Uts, fecache)
end

#########
# Utils #
#########
function _make_uh_from_us(op, us, cache)
  Ut, = cache

  u = EvaluationFunction(Ut[1], us[1])

  dus = ()
  for i in 2:get_order(op)+1
    dus = (dus..., EvaluationFunction(Ut[i], us[i]))
  end

  TransientCellField(u, dus)
end
