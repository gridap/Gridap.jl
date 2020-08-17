module CellDataTests

using Test

@testset "CellFields" begin include("CellFieldsTests.jl") end

@testset "CellQuadratures" begin include("CellQuadraturesTests.jl") end

@testset "QPointCellFields" begin include("QPointCellFieldsTests.jl") end

end # module
