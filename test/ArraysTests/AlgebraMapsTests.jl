module AlgebraMapsTests

using Test
using Gridap.Arrays
using LinearAlgebra

a = rand(3,4)
b = rand(4,5)

c = a*b
test_map(c,*,a,b)
cache = return_cache(*,a,b)
@test evaluate!(cache,*,a,b) === evaluate!(cache,*,a,b)

a = rand(3,4)
b = rand(4)

c = a*b
test_map(c,*,a,b)
cache = return_cache(*,a,b)
@test evaluate!(cache,*,a,b) === evaluate!(cache,*,a,b)

a = rand(3,4)
b = 2.0

c = a*b
test_map(c,*,a,b)
cache = return_cache(*,a,b)
@test evaluate!(cache,*,a,b) === evaluate!(cache,*,a,b)

c = b*a
test_map(c,*,b,a)
cache = return_cache(*,b,a)
@test evaluate!(cache,*,b,a) === evaluate!(cache,*,b,a)

a = 3
b = 4

c = a*b
test_map(c,*,a,b)

k! = AddEntriesMap(+)

A = zeros(4,5)
k!(A,fill(3.0,2,2),[1,3],[2,1])
k!(A,fill(3.0,2,2),[1,3],[2,1])
@test A[[1,3],[2,1]] == fill(6.0,2,2)

A = zeros(4)
k!(A,fill(3.0,2),[1,3])
k!(A,fill(3.0,2),[1,3])
@test A[[1,3]] == fill(6.0,2)

k! = TouchEntriesMap()
A = zeros(4,5)
k!(A,fill(3.0,2,2),[1,1],[2,1])
@test A == zeros(4,5)

A = zeros(4)
k!(A,fill(3.0,2),[1,1])
@test A == zeros(4)

# LocalStagedSolveMap

@testset "LocalStagedSolveMap" begin
  # Blocks: tests (2,3,1), trials (2,4), inputs (3,2). Lower triangular in the stages
  # ((1,),(1,)) and ((2,3),(2,)).
  rows, cols, icols = (2, 3, 1), (2, 4), (3, 2)
  A = [randn(r, c) for r in rows, c in cols]
  A[1, 2] .= 0
  A[1, 1] += 5*I
  B = [randn(r, c) for r in rows, c in icols]
  lhs = ArrayBlock(A, fill(true, 3, 2))
  rhs = ArrayBlock(B, fill(true, 3, 2))
  X = reduce(vcat, [reduce(hcat, A[i, :]) for i in 1:3]) \ reduce(vcat, [reduce(hcat, B[i, :]) for i in 1:3])

  staged = LocalStagedSolveMap((((1,), (1,)), ((2, 3), (2,))))
  mono = LocalStagedSolveMap((((1, 2, 3), (1, 2)),))
  first_block = LocalStagedSolveMap((((1,), (1,)), ((2, 3), (2,))); output = 1)
  both = LocalStagedSolveMap((((1,), (1,)), ((2, 3), (2,))); output = (1, 2))
  for (k, out) in ((staged, 3:6), (mono, 3:6), (first_block, 1:2), (both, 1:6))
    x = evaluate(k, (lhs, rhs))
    @test length(x) == 2
    @test x[1] ≈ X[out, 1:3]
    @test x[2] ≈ X[out, 4:5]
  end

  # A single right-hand side, with an untouched block
  b = ArrayBlock([randn(2), zeros(3), randn(1)], [true, false, true])
  x = evaluate(staged, (lhs, b))
  xref = reduce(vcat, [reduce(hcat, A[i, :]) for i in 1:3]) \ vcat(b.array[1], zeros(3), b.array[3])
  @test length(x) == 1
  @test only(x) ≈ xref[3:6]
end

end # module
