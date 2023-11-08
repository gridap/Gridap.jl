using Pkg
Pkg.add(url="https://github.com/tamaratambyah/Gridap.jl", rev="rungekutta")

using Gridap
using Plots
using Test

n = 32
p = 1
degree = 4*(p+1)
L = 1
dx = L/n

domain = (0.0, L, 0.0, L)
partition = (n, n)
model  = CartesianDiscreteModel(domain,partition, isperiodic=(true,true))
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

println("FE spaces")
V = TestFESpace(model,
                ReferenceFE(raviart_thomas,Float64,p),
                conformity=:Hdiv)
U = TransientTrialFESpace(V)

W = TestFESpace(model,
                ReferenceFE(lagrangian ,Float64,p),
                conformity=:L2)
R = TransientTrialFESpace(W)

X = TransientMultiFieldFESpace([U,R])
Y = MultiFieldFESpace([V,W])


g = 1.0
H0 = 1.0

# initial condition
u0(x,t) = VectorValue(0.0, 0.0)
u0(t) = x -> u0(x,t)
h0(x,t) = H0 + 0.1*cos(2*π*x[1])
h0(t) = x -> h0(x,t)

uh0 = interpolate_everywhere(u0(0.0),U(0.0))
hh0 = interpolate_everywhere(h0(0.0),R(0.0))
xh0 = interpolate_everywhere([uh0,hh0],X(0.0))


#
lhs(t,(u,h),(v,w)) = ∫( u⋅v + h*w )dΩ
rhs(t,(u,h),(v,w)) = ∫( g*h*(∇⋅v)  )dΩ -  ∫( H0*(∇⋅u)*w  )dΩ
jac(t,(u,h),(du,dh),(v,w)) = ∫( g*(dh*(∇⋅v))  )dΩ -  ∫( (H0*(∇⋅du))*w  )dΩ
jac_t(t,(u,h),(dut,dht),(v,w)) = ∫( dut⋅v + dht*w )dΩ

op = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,X,Y)

t0 = 0.0
tF = 1.0
dt = 0.001


ls = LUSolver()
ode_solver = EXRungeKutta(ls,dt,:EX_SSP_3_0_3)

sol_t = solve(ode_solver,op,xh0,t0,tF)


createpvd("my_tests/transient_sol/wave_eq") do pvd
  pvd[0] = createvtk(Ω,"my_tests/transient_sol/wave_eq_0"*".vtu",cellfields=["u"=>uh0, "h"=>hh0])
  for (xh,t) in sol_t
    println(t)
    uh = xh[1]
    hh = xh[2]
    pvd[t] = createvtk(Ω,"my_tests/transient_sol/wave_eq_$t"*".vtu",cellfields=["u"=>uh, "h"=>hh])
  end
end

global uh = 1
for (xh,t) in sol_t
  println(t)
  # global uh = xh
  # hh = xh[2]
end

Uh = Gridap.ODEs.TransientFETools.allocate_trial_space(sol_t.trial)

Xh0 = get_free_dof_values(xh0)
tf = 0.0
evaluate!(Xh0,sol_t.trial,tf)
sol_t.trial

evaluate!(Uh,sol.trial,tf)


get_free_dof_values(uh)

##### INCLUDE CORILOIS
f = 1.0
function vecPerp(Fold)
  Fnew = VectorValue(-copy(Fold[2]),copy(Fold[1]))
  return Fnew
end

lhs(t,(u,h),(v,w)) = ∫( u⋅v + h*w )dΩ
rhs(t,(u,h),(v,w)) = ∫( g*h*(∇⋅v)  )dΩ  - ∫( f*(vecPerp∘(u)⋅v) )dΩ  -  ∫( H0*(∇⋅u)*w  )dΩ
jac(t,(u,h),(du,dh),(v,w)) = ∫( g*(dh*(∇⋅v))  )dΩ  - ∫( f*(vecPerp∘(du)⋅v) )dΩ -  ∫( (H0*(∇⋅du))*w  )dΩ
jac_t(t,(u,h),(dut,dht),(v,w)) = ∫( dut⋅v + dht*w )dΩ

op = TransientEXRungeKuttaFEOperator(lhs,rhs,jac,jac_t,X,Y)

t0 = 0.0
tF = 1.0
dt = 0.001


ls = LUSolver()
ode_solver = EXRungeKutta(ls,dt,:EX_SSP_3_0_3)

sol_t = solve(ode_solver,op,xh0,t0,tF)


createpvd("my_tests/transient_sol/wave_eq_f") do pvd
  pvd[0] = createvtk(Ω,"my_tests/transient_sol/wave_eq_f_0"*".vtu",cellfields=["u"=>uh0, "h"=>hh0])
  for (xh,t) in sol_t
    println(t)
    uh = xh[1]
    hh = xh[2]
    pvd[t] = createvtk(Ω,"my_tests/transient_sol/wave_eq_f_$t"*".vtu",cellfields=["u"=>uh, "h"=>hh])
  end
end


for (xh,t) in sol_t
  println(t)
  uh = xh[1]
  hh = xh[2]
end
