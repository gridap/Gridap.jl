module DifferentialFormsTests

using Test

@testset "CellMetrics" begin include("CellMetricsTests.jl") end

@testset "DifferentialFormCellFields" begin include("DifferentialFormCellFieldsTests.jl") end

@testset "FeecPoisson" begin include("FeecPoissonTests.jl") end

end # module
