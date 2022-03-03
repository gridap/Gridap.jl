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

function allocate_cache(op::ODEOpFromFEOp)
  Ut = get_trial(op.feop)
  U = allocate_trial_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for i in 1:get_order(op)
    Uts = (Uts...,∂t(Uts[i]))
    Us = (Us...,allocate_trial_space(Uts[i+1]))
  end
  fecache = allocate_cache(op.feop)
  ode_cache = (Us,Uts,fecache)
  ode_cache
end

function allocate_cache(op::ODEOpFromFEOp,v::AbstractVector,a::AbstractVector)
  ode_cache = allocate_cache(op)
  (v,a, ode_cache)
end

function update_cache!(ode_cache,op::ODEOpFromFEOp,t::Real)
  _Us,Uts,fecache = ode_cache
  Us = ()
  for i in 1:get_order(op)+1
    Us = (Us...,evaluate!(_Us[i],Uts[i],t))
  end
  fecache = update_cache!(fecache,op.feop,t)
  (Us,Uts,fecache)
end

function allocate_residual(op::ODEOpFromFEOp,uhF::AbstractVector,ode_cache)
  Us,Uts,fecache = ode_cache
  uh = EvaluationFunction(Us[1],uhF)
  allocate_residual(op.feop,uh,fecache)
end

function allocate_jacobian(op::ODEOpFromFEOp,uhF::AbstractVector,ode_cache)
  Us,Uts,fecache = ode_cache
  uh = EvaluationFunction(Us[1],uhF)
  allocate_jacobian(op.feop,uh,fecache)
end

"""
It provides A(t,uh,∂tuh,...,∂t^Nuh) for a given (t,uh,∂tuh,...,∂t^Nuh)
"""
function residual!(
  b::AbstractVector,
  op::ODEOpFromFEOp,
  t::Real,
  xhF::Tuple{Vararg{AbstractVector}},
  ode_cache)
  Xh, = ode_cache
  dxh = ()
  for i in 2:get_order(op)+1
    dxh = (dxh...,EvaluationFunction(Xh[i],xhF[i]))
  end
  xh=TransientCellField(EvaluationFunction(Xh[1],xhF[1]),dxh)
  residual!(b,op.feop,t,xh,ode_cache)
end


"""
It adds contribution to the Jacobian with respect to the i-th time derivative,
with i=0,...,N. That is, adding γ_i*[∂A/∂(∂t^iuh)](t,uh,∂tuh,...,∂t^Nuh) for a
given (t,uh,∂tuh,...,∂t^Nuh) to a given matrix J, where γ_i is a scaling coefficient
provided by the `ODESolver`, e.g., 1/Δt for Backward Euler; It represents
∂(δt^i(uh))/∂(uh), in which δt^i(⋅) is the approximation of ∂t^i(⋅) in the solver.
Note that for i=0, γ_i=1.0.
"""
function jacobian!(
  A::AbstractMatrix,
  op::ODEOpFromFEOp,
  t::Real,
  xhF::Tuple{Vararg{AbstractVector}},
  i::Integer,
  γᵢ::Real,
  ode_cache)
  Xh, = ode_cache
  dxh = ()
  for i in 2:get_order(op)+1
    dxh = (dxh...,EvaluationFunction(Xh[i],xhF[i]))
  end
  xh=TransientCellField(EvaluationFunction(Xh[1],xhF[1]),dxh)
  jacobian!(A,op.feop,t,xh,i,γᵢ,ode_cache)
end

"""
Add the contribution of all jacobians ,i.e., ∑ᵢ γ_i*[∂A/∂(∂t^iuh)](t,uh,∂tuh,...,∂t^Nuh,vh)
"""
function jacobians!(
  J::AbstractMatrix,
  op::ODEOpFromFEOp,
  t::Real,
  xhF::Tuple{Vararg{AbstractVector}},
  γ::Tuple{Vararg{Real}},
  ode_cache)
  Xh, = ode_cache
  dxh = ()
  for i in 2:get_order(op)+1
    dxh = (dxh...,EvaluationFunction(Xh[i],xhF[i]))
  end
  xh=TransientCellField(EvaluationFunction(Xh[1],xhF[1]),dxh)
  jacobians!(J,op.feop,t,xh,γ,ode_cache)
end
