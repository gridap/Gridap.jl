module MultiFEFunctionsTests

using Test

using Gridap

order = 1
diritag = "boundary"
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(2,3))
fespace = H1ConformingFESpace(Float64,model,order,diritag)

ufun1(x) = x[1] + x[2]
U1 = TrialFESpace(fespace,ufun1)

ufun2(x) = x[1] + x[2]
U2 = TrialFESpace(fespace,ufun2)

U = MultiFESpace([U1,U2])

n = num_free_dofs(U)

x = rand(n)

uh = MultiFEFunction(x,U)

@test length(uh) == 2
@test free_dofs(uh[1]) == x[1:2]
@test free_dofs(uh[2]) == x[3:4]

uh1, state = iterate(uh)
@test free_dofs(uh1) == x[1:2]

uh2, state = iterate(uh,state)
@test free_dofs(uh2) == x[3:4]

zh = zero(U)
@test isa(zh,MultiFEFunction)
@test free_dofs(zh) == zeros(num_free_dofs(U))

zh = FEFunction(U,free_dofs(zh))
@test isa(zh,MultiFEFunction)
@test free_dofs(zh) == zeros(num_free_dofs(U))

end # module MultiFEFunctionsTests
