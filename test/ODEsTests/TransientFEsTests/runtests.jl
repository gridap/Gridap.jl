module TransientFEToolsTests

using Test

@testset "TransientFETests" begin include("TransientFETests.jl") end

@testset "TransientFEOperatorsTests" begin include("TransientFEOperatorsTests.jl") end

@testset "Transient2ndOrderFEOperatorsTests" begin include("Transient2ndOrderFEOperatorsTests.jl") end

@testset "AffineFEOperatorsTests" begin include("AffineFEOperatorsTests.jl") end

@testset "ConstantFEOperatorsTests" begin include("ConstantFEOperatorsTests.jl") end

@testset "HeatEquationTests" begin include("HeatEquationTests.jl") end

@testset "HeatVectorEquationTests" begin include("HeatVectorEquationTests.jl") end

@testset "VectorHeatEquationTests" begin include("VectorHeatEquationTests.jl") end

@testset "StokesEquationTests" begin include("StokesEquationTests.jl") end

@testset "BoundaryEquationTests" begin include("BoundaryHeatEquationTests.jl") end

@testset "DGHeatEquationTests" begin include("DGHeatEquationTests.jl") end

@testset "FreeSurfacePotentialFlowTests" begin include("FreeSurfacePotentialFlowTests.jl") end

@testset "HeatEquationAutoDiffTests" begin include("HeatEquationAutoDiffTests.jl") end

@testset "StokesEquationAutoDiffTests" begin include("StokesEquationAutoDiffTests.jl") end

@testset "ForwardEulerHeatEquationTests" begin include("ForwardEulerHeatEquationTests.jl") end

end # module
