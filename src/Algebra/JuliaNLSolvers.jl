module JuliaNLSolvers

using Gridap
using NLsolve

export JuliaNLSolver

import Gridap: solve!

mutable struct JuliaNLSolver <: NonLinearSolver
  ls::LinearSolver
  kwargs::Dict
end

function JuliaNLSolver(ls::LinearSolver;kwargs...)
  @assert ! haskey(kwargs,:linsolve) "linsolve cannot be used here. It is managed internally"
  JuliaNLSolver(ls,kwargs)
end

function JuliaNLSolver(;kwargs...)
  ls = BackslashSolver()
  JuliaNLSolver(ls;kwargs...)
end

mutable struct JuliaNLSolversCache
  df::OnceDifferentiable
  ss::SymbolicSetup
  result
end

function solve!(x::AbstractVector,nls::JuliaNLSolver,op::NonLinearOperator)
  cache = _setup_cache(x,nls,op)
  solve!(x,nls,op,cache)
  cache
end

function solve!(
  x::AbstractVector,nls::JuliaNLSolver,op::NonLinearOperator,cache::JuliaNLSolversCache)
  df = cache.df
  ss = cache.ss
  kwargs = nls.kwargs
  function linsolve!(x,A,b)
    ns = numerical_setup(ss,A)
    solve!(x,ns,A,b)
  end
  r = nlsolve(df,x;linsolve=linsolve!,kwargs...)
  cache.result = r
  x[:] .= r.zero
end

function _setup_cache(x0,nls,op)

  f!(r,x) = residual!(r,op,x)
  j!(j,x) = jacobian!(j,op,x)
  
  f0 = residual(op,x0)
  j0 = jacobian(op,x0)
  
  df = OnceDifferentiable(f!,j!,x0,f0,j0)

  ss = symbolic_setup(nls.ls,j0)

  JuliaNLSolversCache(df,ss,nothing)
end

end # module
