module TensorValuesTests

using Test

@testset "TypesTests" begin
  include("TypesTests.jl")
end

@testset "OperationsTests" begin
  include("OperationsTests.jl")
end

@testset "IndexingTests" begin
  include("IndexingTests.jl")
end

@testset "ReinterpretTests" begin
  include("ReinterpretTests.jl")
end

end # module
