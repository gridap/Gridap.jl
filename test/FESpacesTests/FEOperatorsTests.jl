module FEOperatorsTests

##
using Test
using Gridap

import Gridap: ∇
##

# Define manufactured functions
ufun(x) = x[1] + x[2]
ufun_grad(x) = VectorValue(1.0,1.0)
∇(::typeof(ufun)) = ufun_grad
bfun(x) = 0.0

# Construct the discrete model
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

order = 1

include("FEOperatorsTestsMixin.jl")

include("FEOperatorsTestsRobin.jl")

order = 3

include("FEOperatorsTestsMixin.jl")

include("FEOperatorsTestsRobin.jl")

model = simplexify(model)

order = 1

include("FEOperatorsTestsMixin.jl")

include("FEOperatorsTestsRobin.jl")

order = 3

include("FEOperatorsTestsMixin.jl")

include("FEOperatorsTestsRobin.jl")


# Further tests
op = LinearFEOperator(V,U,t_Ω,t_ΓN,t_ΓR,t_ΓS)
op = LinearFEOperator(V,U,t_Ω)

end # module
