"""

The exported names are
$(EXPORTS)
"""
module DiffEqsWrappers

using DocStringExtensions

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

  odeop = get_algebraic_operator(op)
  odeopcache = allocate_cache(odeop)

  function _residual!(res, du, u, p, t)
    # TO DO (minor): Improve update_cache! st do nothing if same time t as in the cache
    # now it would be done twice (residual and jacobian)
    odeopcache = update_cache!(odeopcache, odeop, t)
    residual!(res, odeop, t, (u, du), odeopcache)
  end

  function _jacobian!(jac, du, u, p, gamma, t)
    odeopcache = update_cache!(odeopcache, odeop, t)
    z = zero(eltype(jac))
    fillstored!(jac, z)
    jacobians!(jac, odeop, t, (u, du), (1.0, gamma), odeopcache)
  end

  function _mass!(mass, du, u, p, t)
    odeopcache = update_cache!(odeopcache, odeop, t)
    z = zero(eltype(mass))
    fillstored!(mass, z)
    jacobian!(mass, odeop, t, (u, du), 2, 1.0, odeopcache)
  end

  function _stiffness!(stif, du, u, p, t)
    odeopcache = update_cache!(odeopcache, odeop, t)
    z = zero(eltype(stif))
    fillstored!(stif, z)
    jacobian!(stif, odeop, t, (u, du), 1, 1.0, odeopcache)
  end

  return _residual!, _jacobian!, _mass!, _stiffness!

end

"""
  It allocates the Jacobian (or mass or stiffness) matrix, given the `FEOperator`
  and a vector of size total number of unknowns
"""
function prototype_jacobian(op::TransientFEOperator, u0)
  odeop = get_algebraic_operator(op)
  odeopcache = allocate_cache(odeop) # Not acceptable in terms of performance
  return allocate_jacobian(odeop, u0, odeopcache)
end

const prototype_mass = prototype_jacobian

const prototype_stiffness = prototype_jacobian

end # module DiffEqsWrappers
