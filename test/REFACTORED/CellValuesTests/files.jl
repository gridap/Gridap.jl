
include("CellValuesMocks.jl")

include("MapsMocks.jl")

@testset "CellValuesMocks" begin include("CellValuesMocksTests.jl") end

@testset "MapsMocks" begin include("MapsMocksTests.jl") end

@testset "CellValuesGallery" begin include("CellValuesGalleryTests.jl") end

@testset "Kernels" begin include("KernelsTests.jl") end

@testset "CellNumberApply" begin include("CellNumberApplyTests.jl") end

@testset "CellArrayApply" begin include("CellArrayApplyTests.jl") end

@testset "MapApply" begin include("MapApplyTests.jl") end

@testset "CellMapApply" begin include("CellMapApplyTests.jl") end

