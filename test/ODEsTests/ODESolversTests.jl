module ODESolversTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")
include("ODESolversMocks.jl")

num_eqs = 5
order_max = 5

all_mats = ntuple(_ -> randn(num_eqs, num_eqs), order_max + 1)
all_forms = ntuple(k -> (t -> all_mats[k] .* cospi(t)), order_max + 1)

vec = randn(num_eqs)
forcing(t) = vec .* cospi(t)

mat0 = zeros(num_eqs, num_eqs)
form0(t) = mat0

vec0 = zeros(num_eqs)
forcing0(t) = vec0

t0 = randn()
tF = t0 + rand()
dt = (tF - t0) / 10
all_us0 = ntuple(i -> randn(num_eqs), order_max + 1)
all_usF = copy.(all_us0)
exp_usF = copy.(all_us0)

r = zeros(num_eqs)
J = spzeros(num_eqs, num_eqs)

atol = 1.0e-12
rtol = 1.0e-8
maxiter = 100
disslvr = DiscreteODESolverMock(rtol, atol, maxiter)
odeslvr = ODESolverMock(disslvr, dt)

for N in 1:order_max
  us0 = tuple((all_us0[k] for k in 1:N)...)
  usF = tuple((all_usF[k] for k in 1:N)...)
  forms = all_forms[1:N+1]

  # Compute the expected solution after one step of the ODE solver
  tx = t0 + dt

  f = forcing(tx)
  fill!(r, zero(eltype(r)))
  r .+= f

  m = last(forms)(tx)
  fillstored!(J, zero(eltype(J)))
  J .+= m

  ws = ntuple(i -> 0, N + 1)
  ws = Base.setindex(ws, 1, N + 1)
  for i in N:-1:1
    ui0, uiF = us0[i], usF[i]
    copy!(uiF, ui0)
    wi, coef = 0, 1
    for j in i+1:N
      coef = coef * dt / (j - i)
      wi += coef * ws[j]
      axpy!(coef, usF[j], uiF)
    end
    coef = coef * dt / (N + 1 - i)
    wi += coef * ws[N+1]
    ws = Base.setindex(ws, wi, i)

    usF = Base.setindex(usF, uiF, i)

    form = forms[i](tx)
    r .+= form * uiF
    J .+= wi .* form
  end

  # Solve system
  rmul!(r, -1)
  x = J \ r

  # Finalize and store solution in exp_usF
  for i in N:-1:1
    global exp_usF
    ui0, exp_uiF = us0[i], exp_usF[i]
    copy!(exp_uiF, ui0)
    coef = 1
    for j in i+1:N
      coef = coef * dt / (j - i)
      axpy!(coef, exp_usF[j], exp_uiF)
    end
    coef = coef * dt / (N + 1 - i)
    axpy!(coef, x, exp_uiF)
    exp_usF = Base.setindex(exp_usF, exp_uiF, i)
  end

  for C in (NonlinearODE, QuasilinearODE, SemilinearODE, LinearODE,)
    standard_odeop = ODEOperatorMock{C}(forms, forcing)

    # Create an IMEXODEOperator randomly
    im_forms = ()
    ex_forms = ()
    for k in 0:N-1
      form = forms[k+1]
      to_im = rand(Bool)
      im_forms = (im_forms..., to_im ? form : form0)
      ex_forms = (ex_forms..., to_im ? form0 : form)
    end
    im_forms = (im_forms..., last(forms))

    to_im = rand(Bool)
    im_forcing = to_im ? forcing : forcing0
    ex_forcing = to_im ? forcing0 : forcing

    im_odeop = ODEOperatorMock{C}(im_forms, im_forcing)
    ex_odeop = ODEOperatorMock{C}(ex_forms, ex_forcing)
    imex_odeop = IMEXODEOperator(im_odeop, ex_odeop)

    for odeop in (standard_odeop, imex_odeop,)
      tF, usF, cache = solve_step!(usF, odeslvr, odeop, t0, us0)

      for i in 1:N
        @test usF[i] ≈ exp_usF[i]
      end

      @test test_ode_solver(odeslvr, odeop, t0, us0, t0, us0, dt, tF)
    end
  end
end

end # module ODESolversTests
