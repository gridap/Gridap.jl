module ODESolversTests

using Test
using SparseArrays

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")
include("ODESolversMocks.jl")

num_eqs = 5
order_max = 5

all_mats = ntuple(_ -> sprandn(num_eqs, num_eqs, 1.0), order_max + 1)
all_forms = ntuple(k -> (t -> all_mats[k] .* cospi(t)), order_max + 1)

vec = randn(num_eqs)
forcing(t) = vec .* cospi(t)

mat0 = sprand(num_eqs, num_eqs, 1.0)
nonzeros(mat0) .= 0
form0(t) = mat0

vec0 = zeros(num_eqs)
forcing0(t) = vec0

t0 = randn()
tF = t0 + rand()
dt = (tF - t0) / 10
all_us0 = ntuple(i -> randn(num_eqs), order_max + 1)
all_usF = copy.(all_us0)
exp_usF = copy.(all_us0)

exp_J = spzeros(num_eqs, num_eqs)
exp_r = zeros(num_eqs)

atol = 1.0e-12
rtol = 1.0e-8
maxiter = 100
sysslvr = NonlinearSolverMock(rtol, atol, maxiter)
odeslvr = ODESolverMock(sysslvr, dt)

for N in 1:order_max
  us0 = ntuple(k -> all_us0[k], N)
  usF = ntuple(k -> all_usF[k], N)
  forms = all_forms[1:N+1]

  # Compute the expected solution after one step of the ODE solver
  tx = t0 + dt

  f = forcing(tx)
  fill!(exp_r, zero(eltype(exp_r)))
  exp_r .= -f

  m = last(forms)(tx)
  fillstored!(exp_J, zero(eltype(exp_J)))
  exp_J .+= m

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
    exp_r .+= form * uiF
    exp_J .+= wi .* form
  end

  # Solve system
  rmul!(exp_r, -1)
  exp_x = exp_J \ exp_r

  # Update state
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
    axpy!(coef, exp_x, exp_uiF)
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
      # Allocate cache
      odecache = allocate_odecache(odeslvr, odeop, t0, us0)

      # Starting procedure
      state0, odecache = ode_start(
        odeslvr, odeop,
        t0, us0,
        odecache
      )

      # Marching procedure
      stateF = copy.(state0)
      tF, stateF, odecache = ode_march!(
        stateF,
        odeslvr, odeop,
        t0, state0,
        odecache
      )

      # Finishing procedure
      uF = copy(first(us0))
      uF, odecache = ode_finish!(
        uF,
        odeslvr, odeop,
        t0, tF, stateF,
        odecache
      )

      usF = stateF
      for i in 1:N
        @test usF[i] ≈ exp_usF[i]
      end
      @test uF ≈ first(exp_usF)

      @test test_ode_solver(odeslvr, odeop, t0, us0)
    end
  end
end

end # module ODESolversTests
