module ODEOperatorsTests

using Test
using SparseArrays

using Gridap
using Gridap.ODEs

include("ODEOperatorsMocks.jl")

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

t = randn()
all_us = ntuple(i -> randn(num_eqs), order_max + 1)

exp_r = zeros(num_eqs)
exp_J = spzeros(num_eqs, num_eqs)

Ts = (NonlinearODE, QuasilinearODE, SemilinearODE, LinearODE)

for N in 0:order_max
  us = tuple((all_us[k] for k in 1:N+1)...)
  forms = all_forms[1:N+1]

  for T in Ts
    standard_odeop = ODEOperatorMock{T}(forms, forcing)
    odeops = (standard_odeop,)

    # Randomly create a `IMEXODEOperator`s
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

      im_odeop = ODEOperatorMock{T}(im_forms, im_forcing)
      for T_ex in Ts
        ex_odeop = ODEOperatorMock{T_ex}(ex_forms, ex_forcing)
        imex_odeop = IMEXODEOperator(im_odeop, ex_odeop)
        odeops = (odeops..., imex_odeop)
      end
    end

    # Compute expected residual
    f = forcing(t)
    copy!(exp_r, -f)
    for (ui, formi) in zip(us, forms)
      form = formi(t)
      exp_r .+= form * ui
    end

    for odeop in odeops
      num_forms = get_num_forms(odeop)
      if odeop isa IMEXODEOperator
        im_odeop, ex_odeop = get_imex_operators(odeop)
        T_im, T_ex = ODEOperatorType(im_odeop), ODEOperatorType(ex_odeop)
        if T_im <: AbstractLinearODE
          if T_ex <: AbstractLinearODE
            @test num_forms == get_order(im_odeop) + 1
          else
            @test num_forms == 1
          end
        elseif T_im <: AbstractQuasilinearODE
          @test num_forms == 1
        else
          @test num_forms == 0
        end
      else
        if T <: AbstractLinearODE
          @test num_forms == get_order(odeop) + 1
        elseif T <: AbstractQuasilinearODE
          @test num_forms == 1
        else
          @test num_forms == 0
        end
      end

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
      end
      ws = ntuple(_ -> 1, N + 1)
      jacobian!(J, odeop, t, us, ws, odeopcache)
      @test J ≈ exp_J

      @test test_ode_operator(odeop, t, us)
    end
  end
end

end # module ODEOperatorsTests
