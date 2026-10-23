module DifferentialFormsTests

using Test

@testset "CellMetrics" begin include("CellMetricsTests.jl") end

@testset "MetricOperations" begin include("MetricOperationsTests.jl") end

@testset "MetricCodifferential" begin include("MetricCodifferentialTests.jl") end

@testset "MeasureConvention" begin include("MeasureConventionTests.jl") end

@testset "MetricFormAssembly" begin include("MetricFormAssemblyTests.jl") end

@testset "MaterialHodgeStars" begin include("MaterialHodgeStarsTests.jl") end

@testset "InclusionTraces" begin include("InclusionTracesTests.jl") end

@testset "DifferentialFormCellFields" begin include("DifferentialFormCellFieldsTests.jl") end

@testset "FeecPoisson" begin include("FeecPoissonTests.jl") end

end # module
