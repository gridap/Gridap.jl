module ODEOperatorsTests

using Test

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")

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

t = randn()
all_us = ntuple(i -> randn(num_eqs), order_max + 1)

exp_r = zeros(num_eqs)
exp_J = spzeros(num_eqs, num_eqs)

for N in 0:order_max
  us = tuple((all_us[k] for k in 1:N+1)...)
  forms = all_forms[1:N+1]

  for C in (NonlinearODE, QuasilinearODE, SemilinearODE, LinearODE)
    standard_odeop = ODEOperatorMock{C}(forms, forcing)
    odeops = (standard_odeop,)

    # Create an IMEXODEOperator randomly
    if N > 0
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

      odeops = (odeops..., imex_odeop)
    end

    # Compute expected residual
    f = forcing(t)
    copy!(exp_r, f)
    for (ui, formi) in zip(us, forms)
      form = formi(t)
      exp_r .+= form * ui
    end

    for odeop in odeops
      odeopcache = allocate_odeopcache(odeop, t, us)
      update_odeopcache!(odeopcache, odeop, t)

      r = allocate_residual(odeop, t, us, odeopcache)
      @test size(r) == (num_eqs,)

      J = allocate_jacobian(odeop, t, us, odeopcache)
      @test size(J) == (num_eqs, num_eqs)

      residual!(r, odeop, t, us, odeopcache)
      @test r ≈ exp_r

      fill!(exp_J, zero(eltype(exp_J)))
      fill!(J, zero(eltype(J)))
      for k in 0:N
        exp_J .+= forms[k+1](t)
        jacobian!(J, odeop, t, us, k, 1, odeopcache)
        @test J ≈ exp_J
      end

      @test test_ode_operator(odeop, t, us)
    end
  end
end

end # module ODEOperatorsTests
