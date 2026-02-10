module InverseTests  

using Gridap, Test, ChainRulesCore

reffe = ReferenceFE(lagrangian, Float64, 2)
model = CartesianDiscreteModel((0,1,0,1),(2,2))
V0 = TestFESpace(model, reffe, dirichlet_tags=[1,2,3,4,5,7,8])
U0 = TrialFESpace(V0, 0.0)

points = [VectorValue(0.4324,0.1232141)]
obs_operator = FEObservationOperator(points,U0)

u = rand(num_free_dofs(U0))

l,pullback =  rrule(obs_operator, u)
du = pullback([1])[2]

function cdf_gradient(f, x, h=1e-4)
  manu_grad = zeros(length(x))
  x⁺, x⁻ = x[:], x[:]
  for i in eachindex(manu_grad)
    x⁺[i] += h
    x⁻[i] -= h

    manu_grad[i] = (f(x⁺) - f(x⁻)) / (2 * h)

    x⁺[i] -= h
    x⁻[i] += h
  end
  manu_grad
end
du_fd = cdf_gradient(u->obs_operator(u)[1],u)
@test du ≈ du_fd

end