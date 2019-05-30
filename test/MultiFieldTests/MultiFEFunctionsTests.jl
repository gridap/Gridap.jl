include("../../src/MultiField/MultiFEFunctions.jl")

module MultiFEFunctionsTests

using Test

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Assemblers
using Gridap.Geometry.Cartesian
using Gridap.CellMaps
using Gridap.MultiCellArrays
using Gridap.CellQuadratures
using Gridap.CellIntegration

using ..MultiFESpaces
using ..MultiAssemblers
using ..MultiFEFunctions

order = 1
diritag = "boundary"
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(2,3))
fespace = ConformingFESpace(Float64,model,order,diritag)

ufun1(x) = x[1] + x[2]
U1 = TrialFESpace(fespace,ufun1)

ufun2(x) = x[1] + x[2]
U2 = TrialFESpace(fespace,ufun2)

V1 = TestFESpace(fespace)
V2 = V1

V = MultiFESpace([V1,V2])
U = MultiFESpace([U1,U2])

assem = MultiSparseMatrixAssembler(V,U)

n = num_free_dofs(U)

x = rand(n)

uh = MultiFEFunction(x,U,assem)

@test length(uh) == 2
@test free_dofs(uh[1]) == x[1:2]
@test free_dofs(uh[2]) == x[3:4]

uh1, state = iterate(uh)
@test free_dofs(uh1) == x[1:2]

uh2, state = iterate(uh,state)
@test free_dofs(uh2) == x[3:4]

end # module MultiFEFunctionsTests
