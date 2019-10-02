module MultiFieldTests

using Test

@testset "MultiCellArrays" begin include("MultiCellArraysTests.jl") end

@testset "MultiCellMaps" begin include("MultiCellMapsTests.jl") end

@testset "MultiFESpaces" begin include("MultiFESpacesTests.jl") end

@testset "MultiFEFunctions" begin include("MultiFEFunctionsTests.jl") end

@testset "MultiFEBases" begin include("MultiFEBasesTests.jl") end

@testset "MultiAssemblers" begin include("MultiAssemblersTests.jl") end

@testset "MultiFEOperators" begin include("MultiFEOperatorsTests.jl") end

@testset "DGMultiFEOperators" begin include("DGMultiFEOperatorsTests.jl") end

@testset "VectorValuedMultiFEOperators" begin include("VectorValuedMultiFEOperatorsTests.jl") end

@testset "DarcyFEOperators" begin include("DarcyFEOperatorsTests.jl") end

@testset "MultiNonLinearFEOperators" begin include("MultiNonLinearFEOperatorsTests.jl") end

@testset "NSCavityTests" begin include("NSCavityTests.jl") end

end # module
