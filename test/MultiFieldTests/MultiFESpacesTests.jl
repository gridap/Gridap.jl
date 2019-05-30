include("../../src/MultiField/MultiFESpaces.jl")

module MultiFESpacesTests

using Test
using Gridap
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.Geometry.Cartesian
using ..MultiFESpaces

model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))

order = 1
diritag = "boundary"
fespace = ConformingFESpace(Float64,model,order,diritag)

ufun1(x) = x[1] + x[2]
U1 = TrialFESpace(fespace,ufun1)

ufun2(x) = x[1] + x[2]
U2 = TrialFESpace(fespace,ufun2)

U = MultiFESpace([U1,U2])

@test length(U) == 2
@test U[1] === U1
@test U[2] === U2

V1, state = iterate(U)
@test V1 === U1

V2, state = iterate(U,state)
@test V2 === U2

@test num_free_dofs(U) == num_free_dofs(U1) + num_free_dofs(U2)

end # module MultiFESpacesTests
