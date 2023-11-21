using Gridap
using Gridap.Adaptivity

let
  domain = (0.0, 1.0, 0.0, 1.0)
  partition = (2,2)
  model = CartesianDiscreteModel(domain, partition) |> simplexify
  num_iterations = 8
  for iteration in 1:num_iterations
    if iteration == 1
      @time model = refine(model)
    else
      @time model = refine(model.model)
    end
  end
end
