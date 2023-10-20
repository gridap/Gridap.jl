using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap
using LaTeXStrings
using Plots
using JLD


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


#### SET UP FOR SOLVER
u0 = get_free_dof_values( interpolate_everywhere(u(0),U(0.0)) )
solver = RungeKutta(ls,ls,0.001,:BE_1_0_1)
op = TransientRungeKuttaFEOperator(Lhs,Rhs,jac,jac_t,U,V)
t0 = 0.0
cache = nothing
uf .= (u0)

# function solve_step!(uf::AbstractVector,
#   solver::RungeKutta,
#   op::ODEOperator,
#   u0::AbstractVector,
#   t0::Real,
#   cache)

  # Unpack variables
  dt = solver.dt
  s = solver.bt.s
  a = solver.bt.a
  b = solver.bt.b
  c = solver.bt.c
  d = solver.bt.d

  import Gridap.ODEs.ODETools: allocate_cache

  # Create cache if not there
  # if cache === nothing
    ode_cache = allocate_cache(op)
    vi = similar(u0)
    ui = [similar(u0)]
    rhs = similar(u0)
    nls_stage_cache = nothing
    nls_update_cache = nothing
  # else
  #   ode_cache, vi, ui, rhs, nls_stage_cache, nls_update_cache = cache
  # end


  import Gridap.ODEs.ODETools: RungeKuttaNonlinearOperator
  import Gridap.ODEs.ODETools: RungeKuttaStageNonlinearOperator
  import Gridap.FESpaces: get_algebraic_operator

  op1 = get_algebraic_operator(op) # I think this happens outside solve_step!

  # Create RKNL stage operator
  nlop_stage = RungeKuttaStageNonlinearOperator(op1,t0,dt,u0,ode_cache,vi,ui,rhs,0,a)

  # Compute intermediate stages
  # for i in 1:s
    i = 1
    # allocate space to store the RHS at i
    if (length(ui) < i)
      push!(ui,similar(u0))
    end

    import Gridap.ODEs.ODETools: update_cache!
    import Gridap.ODEs.ODETools: update!

    # Update time
    ti = t0 + c[i]*dt
    ode_cache = update_cache!(ode_cache,op,ti)
    update!(nlop_stage,ti,ui,i)

    # if(a[i,i]==0)
    #   # Skip stage solve if a_ii=0 => u_i=u_0, f_i = f_0
    #   @. uf = u0
    # else
      # solve at stage i
    import Gridap.Algebra: solve!

      nls_stage_cache = solve!(uf,solver.nls_stage,nlop_stage,nls_stage_cache)
    # end
    # solve!(x::AbstractVector,ls::LinearSolver,A::AbstractMatrix,b::AbstractVector)
    # Update stage unknown
    @. nlop_stage.ui[i] = uf

  # end

  # Update final time
  tf = t0+dt

  # update final solution
  uf .= u0
  for i in 1:s
    uf = uf + dt*b[i]*nlop_stage.ui[i]
  end

  # # Skip final update if not necessary
  # if !(c[s]==1.0 && a[s,:] == b)

  #   # Create RKNL final update operator
  #   ode_cache = update_cache!(ode_cache,op,tf)
  #   nlop_update = RungeKuttaUpdateNonlinearOperator(op,tf,dt,u0,ode_cache,vi,ui,rhs,s,b)

  #   # solve at final update
  #   nls_update_cache = solve!(uf,solver.nls_update,nlop_update,nls_update_cache)

  # end

  # Update final cache
  cache = (ode_cache, vi, ui, rhs, nls_stage_cache, nls_update_cache)

  return (uf,tf,cache)

# end
