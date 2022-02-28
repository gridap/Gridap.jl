"""

The exported names are
$(EXPORTS)
"""
module DiffEqWrappers

using Test

using Gridap.GridapODEs.TransientFETools: TransientFEOperator

using Gridap.GridapODEs.ODETools: allocate_cache
using Gridap.GridapODEs.ODETools: update_cache!
using Gridap.GridapODEs.ODETools: residual!
using Gridap.GridapODEs.ODETools: jacobians!
using Gridap.GridapODEs.ODETools: jacobian!

using Gridap.Algebra: allocate_jacobian

using Gridap.FESpaces: get_algebraic_operator
using LinearAlgebra: fillstored!

export prototype_jacobian
export prototype_mass
export prototype_stiffness

export diffeq_wrappers

"""
  This method takes a `FEOperator` and returns some methods that can be used
  in `DifferentialEquations.jl`. Assuming we want to solve the nonlinear ODE
  `res(t,u,du) = 0`, we return:
  1. `residual!(res, du, u, p, t)`: It returns the residual (in `res`) at
     `(u,du,t)`, following the signature in `DifferentialEquations`.
     For the moment, we do not support parameters.
  2. `jacobian!(jac, du, u, p, gamma, t)`: Idem for the Jacobian. It returns
     `∂res/∂du*γ + ∂res/∂u`
  3. `mass!(mass, du, u, p, t)`: Idem for the mass matrix. It returns
     `∂res/∂du`
  4. `stiffness!(stif, du, u, p, t)`: Idem for the mass matrix. It returns
     `∂res/∂u`
"""
function diffeq_wrappers(op)

  ode_op = get_algebraic_operator(op)
  ode_cache = allocate_cache(ode_op)

  function _residual!(res, du, u, p, t)
    # TO DO (minor): Improve update_cache! st do nothing if same time t as in the cache
    # now it would be done twice (residual and jacobian)
    ode_cache = update_cache!(ode_cache, ode_op, t)
    residual!(res, ode_op, t, (u, du), ode_cache)
  end

  function _jacobian!(jac, du, u, p, gamma, t)
    ode_cache = update_cache!(ode_cache, ode_op, t)
    z = zero(eltype(jac))
    fillstored!(jac, z)
    jacobians!(jac, ode_op, t, (u, du), (1.0, gamma), ode_cache)
  end

  function _mass!(mass, du, u, p, t)
    ode_cache = update_cache!(ode_cache, ode_op, t)
    z = zero(eltype(mass))
    fillstored!(mass, z)
    jacobian!(mass, ode_op, t, (u, du), 2, 1.0, ode_cache)
  end

  function _stiffness!(stif, du, u, p, t)
    ode_cache = update_cache!(ode_cache, ode_op, t)
    z = zero(eltype(stif))
    fillstored!(stif, z)
    jacobian!(stif, ode_op, t, (u, du), 1, 1.0, ode_cache)
  end

  return _residual!, _jacobian!, _mass!, _stiffness!

end

"""
  It allocates the Jacobian (or mass or stiffness) matrix, given the `FEOperator`
  and a vector of size total number of unknowns
"""
function prototype_jacobian(op::TransientFEOperator,u0)
  ode_op = get_algebraic_operator(op)
  ode_cache = allocate_cache(ode_op) # Not acceptable in terms of performance
  return allocate_jacobian(ode_op, u0, ode_cache)
end

const prototype_mass = prototype_jacobian

const prototype_stiffness = prototype_jacobian

end #module
