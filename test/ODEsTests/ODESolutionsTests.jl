module ODESolutionsTests

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
dt = (tF - t0) / 5
all_us0 = ntuple(i -> randn(num_eqs), order_max + 1)
all_usF = copy.(all_us0)
exp_usF = copy.(all_us0)

atol = 1.0e-12
rtol = 1.0e-8
maxiter = 100
sysslvr = NonlinearSolverMock(rtol, atol, maxiter)
odeslvr = ODESolverMock(sysslvr, dt)

Ts = (NonlinearODE, QuasilinearODE, SemilinearODE, LinearODE)

for N in 1:order_max
  us0 = tuple((all_us0[k] for k in 1:N)...)
  usF = tuple((all_usF[k] for k in 1:N)...)
  forms = all_forms[1:N+1]

  for T in Ts
    standard_odeop = ODEOperatorMock{T}(forms, forcing)
    odeops = (standard_odeop,)

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

    for T_ex in Ts
      imex_odeop = IMEXODEOperator(
        ODEOperatorMock{T}(im_forms, im_forcing),
        ODEOperatorMock{T_ex}(ex_forms, ex_forcing)
      )
      odeops = (odeops..., imex_odeop)
    end

    for odeop in odeops
      odesltn = solve(odeslvr, odeop, t0, tF, us0)

      tprev = t0
      for (t_n, u_n) in odesltn
        @test t_n ≈ tprev + dt
        tprev = t_n
      end

      @test test_ode_solution(odesltn)
    end
  end
end

end # module ODESolutionsTests
