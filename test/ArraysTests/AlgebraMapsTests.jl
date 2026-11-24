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

# LocalConstrainedSolveMap

@testset "LocalConstrainedSolveMap" begin
  np, nλ, icols = 4, 2, (3, 2)
  App = randn(np, np) + 5*I
  Apλ, Aλp = randn(np, nλ), randn(nλ, np)
  K = [App Apλ; Aλp zeros(nλ, nλ)]
  lhs = ArrayBlock(reshape([App, Aλp, Apλ, zeros(0, 0)], 2, 2), [true true; true false])
  k = LocalConstrainedSolveMap()

  # Several right-hand sides, one of them without constraint data
  B = reshape([randn(np, icols[1]), randn(nλ, icols[1]), randn(np, icols[2]), zeros(0, 0)], 2, 2)
  rhs = ArrayBlock(B, [true true; true false])
  x = evaluate(k, (lhs, rhs))
  @test isa(x, MatrixBlock) && x.touched == [true true]
  @test x.array[1, 1] ≈ (K \ vcat(B[1, 1], B[2, 1]))[1:np, :]
  @test x.array[1, 2] ≈ (K \ vcat(B[1, 2], zeros(nλ, icols[2])))[1:np, :]

  # A single right-hand side
  b = ArrayBlock([randn(np), randn(nλ)], [true, true])
  x = evaluate(k, (lhs, b))
  @test isa(x, VectorBlock) && x.touched == [true]
  @test x.array[1] ≈ (K \ vcat(b.array[1], b.array[2]))[1:np]

  # Unconstrained 1x1 block system
  lhs1 = ArrayBlock(reshape([App], 1, 1), fill(true, 1, 1))
  rhs1 = ArrayBlock(reshape([B[1, 1], B[1, 2]], 1, 2), fill(true, 1, 2))
  x = evaluate(k, (lhs1, rhs1))
  @test x.touched == [true true]
  @test x.array[1, 1] ≈ App \ B[1, 1]
  @test x.array[1, 2] ≈ App \ B[1, 2]
end

# LocalPenaltySolveMap

@testset "LocalPenaltySolveMap" begin
  # A consistent constraint: the penalised solution is the primal solution of the saddle point
  np, nλ, icols = 4, 1, (3, 2)
  App = Matrix(1.0I, np, np)
  Apλ = randn(np, nλ)
  Aλp = Matrix(Apλ')
  K = [App Apλ; Aλp zeros(nλ, nλ)]
  blocks() = ArrayBlock(reshape([copy(App), copy(Aλp), copy(Apλ), zeros(0, 0)], 2, 2), [true true; true false])
  k = LocalPenaltySolveMap()

  # The inputs are overwritten
  Bp = [Apλ * randn(nλ, n) for n in icols]
  rhs = ArrayBlock(reshape([copy(Bp[1]), zeros(nλ, icols[1]), copy(Bp[2]), zeros(0, 0)], 2, 2), [true true; true false])
  x = evaluate(k, (blocks(), rhs))
  @test isa(x, MatrixBlock) && x.touched == [true true]
  @test size(x.array[1, 1]) == (np, icols[1]) && size(x.array[1, 2]) == (np, icols[2])

  b = ArrayBlock([Apλ * randn(nλ), zeros(nλ)], [true, true])
  x = evaluate(k, (blocks(), ArrayBlock(copy.(b.array), copy(b.touched))))
  @test isa(x, VectorBlock) && x.touched == [true]
  @test length(x.array[1]) == np
end

# LocalSolveMap

@testset "LocalSolveMap" begin
  k = LocalSolveMap()
  A, B = randn(4, 4) + 5*I, randn(4, 3)
  @test evaluate(k, (A, B)) ≈ A \ B

  # Block systems: tests (2,3,1), trials (2,4), inputs (3,2)
  rows, cols, icols = (2, 3, 1), (2, 4), (3, 2)
  A = [randn(r, c) for r in rows, c in cols]
  A[1, 1] += 5*I
  B = [randn(r, c) for r in rows, c in icols]
  lhs = ArrayBlock(A, fill(true, 3, 2))
  rhs = ArrayBlock(B, fill(true, 3, 2))
  K = reduce(vcat, [reduce(hcat, A[i, :]) for i in 1:3])
  X = K \ reduce(vcat, [reduce(hcat, B[i, :]) for i in 1:3])

  x = evaluate(k, (lhs, rhs))
  @test isa(x, MatrixBlock) && all(x.touched) && size(x) == (2, 2)
  @test x.array[1, 1] ≈ X[1:2, 1:3] && x.array[1, 2] ≈ X[1:2, 4:5]
  @test x.array[2, 1] ≈ X[3:6, 1:3] && x.array[2, 2] ≈ X[3:6, 4:5]

  # A single right-hand side, with an untouched block
  b = ArrayBlock([randn(2), zeros(3), randn(1)], [true, false, true])
  x = evaluate(k, (lhs, b))
  xref = K \ vcat(b.array[1], zeros(3), b.array[3])
  @test isa(x, VectorBlock) && all(x.touched) && size(x) == (2,)
  @test x.array[1] ≈ xref[1:2] && x.array[2] ≈ xref[3:6]

  # The result is cached
  cache = return_cache(k, lhs, rhs)
  @test evaluate!(cache, k, lhs, rhs) === evaluate!(cache, k, lhs, rhs)
end

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
  # The solution blocks (output blocks x input blocks)
  for (k, out) in ((staged, 3:6), (mono, 3:6), (first_block, 1:2))
    x = evaluate(k, (lhs, rhs))
    @test isa(x, MatrixBlock) && size(x) == (1, 2) && all(x.touched)
    @test x.array[1, 1] ≈ X[out, 1:3]
    @test x.array[1, 2] ≈ X[out, 4:5]
  end

  # Several output blocks
  x = evaluate(both, (lhs, rhs))
  @test size(x) == (2, 2) && all(x.touched)
  @test x.array[1, 1] ≈ X[1:2, 1:3] && x.array[1, 2] ≈ X[1:2, 4:5]
  @test x.array[2, 1] ≈ X[3:6, 1:3] && x.array[2, 2] ≈ X[3:6, 4:5]

  # A single right-hand side, with an untouched block
  b = ArrayBlock([randn(2), zeros(3), randn(1)], [true, false, true])
  x = evaluate(staged, (lhs, b))
  xref = reduce(vcat, [reduce(hcat, A[i, :]) for i in 1:3]) \ vcat(b.array[1], zeros(3), b.array[3])
  @test isa(x, VectorBlock) && x.touched == [true]
  @test x.array[1] ≈ xref[3:6]
  x = evaluate(both, (lhs, b))
  @test isa(x, VectorBlock) && all(x.touched)
  @test x.array[1] ≈ xref[1:2] && x.array[2] ≈ xref[3:6]
end

end # module
