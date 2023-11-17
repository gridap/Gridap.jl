module ODESolversTests

using Gridap.ODEs
using Gridap.ODEs.ODETools: GenericODESolution
using Gridap.ODEs.ODETools: BackwardEuler
using Gridap.ODEs.ODETools: RungeKutta
using Gridap.ODEs.ODETools: IMEXRungeKutta
using Gridap.ODEs.ODETools: EXRungeKutta
using Gridap.ODEs.ODETools: ForwardEuler
using Gridap.ODEs.ODETools: ThetaMethodNonlinearOperator
using Gridap.ODEs.ODETools: GeneralizedAlpha
using Gridap.ODEs.ODETools: solve!
using Gridap.ODEs
using Gridap.ODEs.ODETools
using Gridap
using Test

# using Gridap.Algebra: residual, jacobian

include("ODEOperatorMocks.jl")

op = ODEOperatorMock{Float64,Constant}(1.0,0.0,1.0,1)

include("ODESolverMocks.jl")

t0 = 0.0
tf = 1.0
dt = 0.1
u0 = ones(2)*2

# NonlinearOperator tests

sop = OperatorMock(op,tf,dt,u0)
isa(sop,NonlinearOperator)

ode_cache = allocate_cache(op)

x = zero_initial_guess(sop)
x .+= 1.0
isa(sop,OperatorMock)
isa(x,AbstractVector)
r = allocate_residual(sop,x)
J = allocate_jacobian(sop,x)
residual!(r,sop,x)
jacobian!(J,sop,x)
@test all(r .== [ -11.0 -11.0])
@test all(J .== [ 9.0 0.0; 0.0 9.0])
_r = residual(sop,x)
_J = jacobian(sop,x)
@test all(_r .== [ -11.0 -11.0])
@test all(_J .== [ 9.0 0.0; 0.0 9.0])

# NLSolver tests

nls = NLSolverMock()
cache = solve!(x,nls,sop)
r, J, dx = cache
@test all(r.==_r)
@test all(J.==_J)
@test all(dx.≈11/9)
@test all(x.≈1+11/9)

#ODESolver tests

odesol = ODESolverMock(nls,dt)
uf = copy(u0)
uf.=1.0

uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,nothing)
uf
@test tf==t0+dt
@test all(uf.≈x)

# ODESolutions

tF = 10.0
sol = GenericODESolution(odesol,op,u0,t0,tF)
current, state = Base.iterate(sol)
uf, tf = current
@test tf==t0+dt
@test all(uf.≈x)

# BackwardEulerNonlinearOperator tests

tf = t0+dt
vf = copy(u0)
sop = ThetaMethodNonlinearOperator(op,tf,dt,u0,ode_cache,vf) # See below
x = zero_initial_guess(sop)
x .+= 1.0
r = allocate_residual(sop,x)
J = allocate_jacobian(sop,x)
residual!(r,sop,x)
jacobian!(J,sop,x)
@test all(r .== [ -11.0 -11.0])
@test all(J .== [ 9.0 0.0; 0.0 9.0])
_r = residual(sop,x)
_J = jacobian(sop,x)
@test all(_r .== [ -11.0 -11.0])
@test all(_J .== [ 9.0 0.0; 0.0 9.0])

ls = LUSolver()

# BackwardEuler tests
dt = 0.01
odesol = BackwardEuler(ls,dt)
uf = copy(u0)
uf.=1.0
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# Affine and nonlinear solvers
op = ODEOperatorMock{Float64,Nonlinear}(1.0,0.0,1.0,1)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

op = ODEOperatorMock{Float64,Affine}(1.0,0.0,1.0,1)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# ThetaMethod tests
odesolθ = ThetaMethod(ls,dt,0.5)
ufθ = copy(u0)
ufθ.=1.0
ufθ, tf, cache = solve_step!(ufθ,odesolθ,op,u0,t0,nothing)
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# RK tests
# RK: BE equivalent
# u1-u0 = dt*u1 => u1 = u0/(1-dt) = 2.2222222222222223
# uf-u0 = dt*u1 => uf = u1
odesol = RungeKutta(ls,ls,dt,:BE_1_0_1)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# RK: CN 2nd order
# k1 = u0
# k2 = u0 + dt * 0.5 * u0 + dt * 0.5 * k2
# k2 = u0 * (1+dt*0.5)/(1-dt*0.5)
# un+1 = u0 + dt * 0.5 * u0 + dt * 0.5 * u0 * (1+dt*0.5)/(1-dt*0.5)
# un+1 = u0 * (1+ dt * 0.5 + dt * 0.5* (1+dt*0.5)/(1-dt*0.5))
odesol = RungeKutta(ls,ls,dt,:CN_2_0_2)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# RK: SDIRK 2nd order
# k1 = u0 + dt * 0.25 * k1
# k1 = u0 * 1/(1-dt*0.25)
# k2 = u0 + dt * 0.5 * k1 + dt * 0.25 * k2
# k2 = u0 * 1/(1-dt*0.25) * (1 + dt*0.5/(1-dt*0.25))
# un+1 = u0 + dt * 0.5 * k1 + dt * 0.5 * k2
# un+1 = u0 + dt * 0.5 * u0 * 1/(1-dt*0.25) + dt * 0.5 * u0 * 1/(1-dt*0.25) * (1 + dt*0.5/(1-dt*0.25))
# un+1 = u0 * (1 + dt*0.5/(1-dt*0.25) + dt*0.5/(1-dt*0.25) * (1 + dt*0.5/(1-dt*0.25))
odesol = RungeKutta(ls,ls,dt,:SDIRK_2_0_2)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# RK: TRBDF (2nd order with some 0 on the diagonal)
odesol = RungeKutta(ls,ls,dt,:TRBDF2_3_2_3)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# IMEX RK tests (explicit part = 0)
odesol = IMEXRungeKutta(ls,ls,dt,:IMEX_FE_BE_2_0_1)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# EX-RK: FE equivalent
odesol = EXRungeKutta(ls,dt,:EX_FE_1_0_1)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# EX-RK: SSP equivalent
odesol = EXRungeKutta(ls,dt,:EX_SSP_3_0_3)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# Forward Euler test
odesol = ForwardEuler(ls,dt)
cache = nothing
uf, tf, cache = solve_step!(uf,odesol,op,u0,t0,cache)
@test tf==t0+dt
@test all( (uf.- u0*exp(dt))  .< 1e-3 )
@test test_ode_solver(odesol,op,u0,t0,tf)

# Newmark test
op_const = ODEOperatorMock{Float64,Constant}(1.0,0.0,0.0,2)
op_const_mat = ODEOperatorMock{Float64,ConstantMatrix}(1.0,0.0,0.0,2)
op_affine = ODEOperatorMock{Float64,Affine}(1.0,0.0,0.0,2)
op_nonlinear = ODEOperatorMock{Float64,Nonlinear}(1.0,0.0,0.0,2)
ops = [op_const, op_const_mat, op_affine, op_nonlinear]
ls = LUSolver()
γ = 0.5
β = 0.25
odesol = Newmark(ls,dt,γ,β)
v0 = ones(2)*(β*dt)
a0 = 0.0*ones(2)
for op in ops
  _uf = copy(u0)
  _uf.=1.0
  _vf = copy(v0)
  _af = copy(a0)
  _cache = nothing
  (_uf, _vf, _af), _tf, _cache = solve_step!((_uf,_vf,_af),odesol,op,(u0,v0,a0),t0,_cache)
  aᵧ = γ*_af .+ (1-γ)*a0
  aᵦ = 2*β*_af .+ (1-2*β)*a0
  @test _tf==t0+dt
  @test all(_vf .≈ (v0 + dt*aᵧ))
  @test all(_uf .≈ (u0 + dt*v0 + 0.5*dt^2*aᵦ))
end

# GeneralizedAlpha test
op = ODEOperatorMock{Float64,Nonlinear}(1.0,0.0,1.0,1)
ρ∞ = 1.0 # Equivalent to θ-method with θ=0.5
αf = 1.0/(1.0 + ρ∞)
αm = 0.5 * (3-ρ∞) / (1+ρ∞)
γ = 0.5 + αm - αf
odesolα = GeneralizedAlpha(ls,dt,ρ∞)
ufα = copy(u0)
v0 = 0.0*ones(2)
vf = copy(v0)
ufα.=1.0
(ufα, vf), tf, cache = solve_step!((ufα,vf),odesolα,op,(u0,v0),t0,nothing)
@test tf==t0+dt
@test all(ufα.≈ufθ)
@test all(vf.≈ 1/(γ*dt) * (ufα-u0) + (1-1/γ)*v0)

# GeneralizedAlpha ∂tt test
op = ODEOperatorMock{Float64,Nonlinear}(0.0,0.0,0.0,2)
γ = 0.5
β = 0.25
ρ∞ = 1.0 # Equivalent to Newmark(0.5, 0.25)
odesolN = Newmark(ls,dt,γ,β)
odesolα = GeneralizedAlpha(ls, dt, ρ∞)
u0 = ones(2)*2
v0 = 0.0*ones(2)
a0 = 0.0*ones(2)
ufN = copy(u0)
ufN .= 1.0
vfN = copy(v0)
afN = copy(a0)
(ufN, vfN, afN), tfN, cache =
    solve_step!((ufN,vfN,afN),odesolN,op,(u0,v0,a0),t0,nothing)
u0 = ones(2)*2
v0 = 0.0*ones(2)
a0 = 0.0*ones(2)
ufα = copy(u0)
ufα .= 1.0
vfα = copy(v0)
afα = copy(a0)
(ufα, vfα, afα), tfα, cache =
    solve_step!((ufα,vfα,afα),odesolα,op,(u0,v0,a0),t0,nothing)
@test tfα==tfN
@test sqrt(sum(abs2.(ufα - ufN))) < 1.0e-10
@test sqrt(sum(abs2.(vfα - vfN))) < 1.0e-10
@test sqrt(sum(abs2.(afα - afN))) < 1.0e-10

end #module
