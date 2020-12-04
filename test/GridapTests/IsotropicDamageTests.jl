module IsotropicDamageTests

using Gridap
using LinearAlgebra
using Test

# max imposed displacement
const udx_max = 0.05

const L = 1

# Elastic model
const E = 2.1e4 # Pa
const ν = 0.3 # dim-less
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))
σe(ε) = λ*tr(ε)*one(ε) + 2*μ*ε # Pa
τ(ε) = sqrt(inner(ε,σe(ε))) # Pa^(1/2)

# Damage model
const σ_u = 5.5 # Pa
const r_0 = σ_u / sqrt(E) # Pa^(1/2)
const H = 0.1 # dim-less

function d(r)
  1 - q(r)/r
end

function q(r)
  r_0 + H*(r-r_0)
end

function new_state(r_in,d_in,ε_in)
  τ_in = τ(ε_in)
  if τ_in <= r_in
    r_out = r_in
    d_out = d_in
    damaged = false
  else
    r_out = τ_in
    d_out = d(r_out)
    damaged = true
  end
  damaged, r_out, d_out
end

function σ(ε_in,r_in,d_in)
  _, _, d_out = new_state(r_in,d_in,ε_in)
  (1-d_out)*σe(ε_in)
end

function dσ(dε_in,ε_in,state)
  damaged, r_out, d_out = state
  if ! damaged
    return (1-d_out)*σe(dε_in)
  else
    c_inc = ((q(r_out) - H*r_out)*inner(σe(ε_in),dε_in))/(r_out^3)
    return (1-d_out)*σe(dε_in) - c_inc*σe(ε_in)
  end
end

u(x) = VectorValue( udx_max*x[1]/L, -ν*udx_max*x[2]/L, -ν*udx_max*x[3]/L )

function main(;n,nsteps)

  domain = (0,L,0,L,0,L)
  partition = (n,n,n)
  model = CartesianDiscreteModel(domain,partition)

  labeling = get_face_labeling(model)
  add_tag_from_tags!(labeling,"ux0",[5,7,13,15,17,19,25])
  add_tag_from_tags!(labeling,"ux1",[2,4,6,8,14,16,18,20,26])
  add_tag_from_tags!(labeling,"uxyz0",[1])
  add_tag_from_tags!(labeling,"uxz0",[3])

  order = 1

  V = TestFESpace(
    model,
    ReferenceFE(:Lagrangian,VectorValue{3,Float64},order),
    labels = labeling,
    dirichlet_tags=["ux0","ux1","uxyz0","uxz0"],
    dirichlet_masks=[(true,false,false),(true,false,false),(true,true,true),(true,false,true)])

  degree = 2*order
  Ω = Triangulation(model)
  dΩ = LebesgueMeasure(Ω,degree)

  r = CellState(r_0,dΩ)
  d = CellState(0.0,dΩ)

  nls = NLSolver(show_trace=false, method=:newton)
  solver = FESolver(nls)

  function step(uh_in,factor)

    udx = factor*udx_max
    u0 = VectorValue(0.0,0.0,0.0)
    ud = VectorValue(udx,0.0,0.0)
    U = TrialFESpace(V,[u0,ud,u0,u0])

    res(u,v) = ∫( inner( ε(v), σ∘(ε(u),r,d) ) )*dΩ
    jac(u,du,v) = ∫( inner( ε(v), dσ∘(ε(du),ε(u),new_state∘(r,d,ε(u))) ) )*dΩ

    op = FEOperator(res,jac,U,V)

    free_values = get_free_values(uh_in)
    uh0 = FEFunction(U,free_values)

    uh_out, = solve!(uh0,solver,op)

    update_state!(new_state,r,d,ε(uh_out))

    uh_out
  end

  factors = collect(1:nsteps)*(1/nsteps)
  uh = zero(V)

  for (istep,factor) in enumerate(factors)
    uh = step(uh,factor)
  end

  e = uh - u

  e_l2 = sqrt(sum(∫(e⋅e)*dΩ))
  @test e_l2 < 1.0e-9

end

main(n=5,nsteps=10)


end # module
