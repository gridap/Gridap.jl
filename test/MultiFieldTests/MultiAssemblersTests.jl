include("../../src/MultiField/MultiAssemblers.jl")

module MultiAssemblersTests

using Test

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Assemblers
using Gridap.Geometry.Cartesian

using ..MultiAssemblers

order = 1
diritag = "boundary"
model = CartesianDiscreteModel(domain=(0.0,1.0,0.0,1.0), partition=(4,4))
fespace = ConformingFESpace(Float64,model,order,diritag)

ufun1(x) = x[1] + x[2]
U1 = TrialFESpace(fespace,ufun1)

ufun2(x) = x[1] + x[2]
U2 = TrialFESpace(fespace,ufun2)

V1 = TestFESpace(fespace)
V2 = V1

V = [V1,V2]
U = [U1,U2]

assem = SparseMatrixAssembler(V,U)

@test isa(assem,MultiSparseMatrixAssembler)

end # module MultiAssemblersTests
