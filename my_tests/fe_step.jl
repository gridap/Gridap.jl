using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap



u(x,t) = x[1]*(1-x[1])*t
u(t) = x -> u(x,t)
∂tu = ∂t(u)
f(t) = x -> ∂t(u)(x,t)-Δ(u(t))(x)

n = 4
p = 2
degree = 4*(p+1)
L = 1
dx = L/n
dt = 0.0001
t0 = 0.0
T = 1.0


domain = (0.0, L)
partition = (n)
model  = CartesianDiscreteModel(domain,partition )
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

V = TestFESpace(model,
                ReferenceFE(lagrangian,Float64,p),
                conformity=:H1,
                dirichlet_tags="boundary")
g(x,t::Real) = 0.0
g(t::Real) = x -> g(x,t)
U = TransientTrialFESpace(V,g)



ls = LUSolver()

A(t,u,v) = ∫(( ∇(v)⊙∇(u) ))dΩ
m(u,v) = ∫(v*u)dΩ
B(v,t) = ∫(v*f(t))dΩ
Lhs(t,u,v) = m(u,v)
Rhs(t,u,v) = B(v,t) - A(t,u,v)
res(t,u,v) = A(t,u,v) + m(∂t(u),v) - B(v,t)
jac(t,u,du,v) = A(t,du,v)
jac_t(t,u,dut,v) = m(dut,v)

fe = ForwardEuler(ls,dt)
op_fe = TransientFEOperator(res,jac,jac_t,U,V)
u0 = interpolate_everywhere(u(0),U(0.0))
sol_fe = solve(fe,op_fe,u0,t0,T)

#### SET UP FOR SOLVER
u0 = get_free_dof_values( interpolate_everywhere(u(0),U(0.0)) )
solver = ForwardEuler(ls,0.001)
op = TransientFEOperator(res,jac,jac_t,U,V)
t0 = 0.0
cache = nothing
uf = u0

# if cache === nothing
import Gridap.ODEs.ODETools: allocate_cache
import Gridap.ODEs.ODETools: update_cache!
import Gridap.ODEs.ODETools: update!

  ode_cache = allocate_cache(op)
  vf = similar(u0)
  nl_cache = nothing
# else
  # ode_cache, vf, nl_cache = cache
# end

dt = solver.dt
tf = t0+dt
# The space should have the boundary conditions at tf
ode_cache = update_cache!(ode_cache,op,t0)

import Gridap.ODEs.ODETools: ForwardEulerNonlinearOperator
import Gridap.FESpaces: get_algebraic_operator

op1 = get_algebraic_operator(op)
nlop = ForwardEulerNonlinearOperator(op1,t0,dt,u0,ode_cache,vf)

import Gridap.ODEs.TransientFETools: allocate_residual
import Gridap.ODEs.TransientFETools: allocate_residual
import Gridap.ODEs.ODETools: allocate_residual

nl_cache = solve!(uf,solver.nls,nlop,nl_cache)

cache = (ode_cache, vf, nl_cache)

return (uf,tf,cache)
