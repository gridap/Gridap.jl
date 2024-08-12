module HelpersTests

using Test

@testset "Preferences" begin include("PreferencesTests.jl") end

@testset "Macros" begin include("MacrosTests.jl") end

@testset "HelperFunctions" begin include("HelperFunctionsTests.jl") end

@testset "GridapTypes" begin include("GridapTypesTests.jl") end

end # module
