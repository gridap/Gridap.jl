module DifferentialFormsTests

using Test

@testset "DifferentialFormCellFields" begin include("DifferentialFormCellFieldsTests.jl") end

@testset "FeecPoisson" begin include("FeecPoissonTests.jl") end

end # module
