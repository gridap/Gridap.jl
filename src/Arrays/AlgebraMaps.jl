
# Make Base.:* behave like a Map

return_value(::typeof(*),a::Number,b::AbstractArray{<:Number}) = a*b

function return_cache(::typeof(*),a::Number,b::AbstractArray{<:Number})
  T = typeof(a*testitem(b))
  N = ndims(b)
  CachedArray(T,N)
end

function evaluate!(cache,::typeof(*),a::Number,b::AbstractArray{<:Number})
  setsize!(cache,size(b))
  c = cache.array
  mul!(c,b,a)
  return c
end

return_value(::typeof(*),a::AbstractArray{<:Number},b::Number) = a*b
return_cache(::typeof(*),a::AbstractArray{<:Number},b::Number) = return_cache(*,b,a)
evaluate!(cache,::typeof(*),a::AbstractArray{<:Number},b::Number) = evaluate!(cache,*,b,a)

function return_cache(::typeof(*),a::AbstractArray{<:Number},b::AbstractArray{<:Number})
  T = typeof(testitem(a)*testitem(b))
  N = ifelse(isa(b,AbstractVector),1,2)
  CachedArray(T,N)
end

function evaluate!(cache,::typeof(*),a::AbstractArray{<:Number},b::AbstractArray{<:Number})
  setsize_op!(*,cache,a,b)
  c = cache.array
  mul!(c,a,b)
  return c
end

function return_cache(::typeof(*), a::Diagonal{<:Number}, b::Diagonal{<:Number})
  T = typeof(testitem(a)*testitem(b))
  Diagonal{T}(undef,0)
end

function evaluate!(cache,::typeof(*), a::Diagonal{<:Number}, b::Diagonal{<:Number})
  setsize_op!(*,cache,a,b)
  mul!(cache,a,b)
  return cache
end

function return_value(::typeof(*),a::ArrayBlock{A,2},b::ArrayBlock{B,1}) where {A,B}
  ri = return_value(*,testvalue(A),testvalue(B))
  array = Vector{typeof(ri)}(undef,size(a.array,1))
  touched = fill(false,size(a.array,1))
  ArrayBlock(array,touched)
end

function return_cache(::typeof(*),a::ArrayBlock,b::ArrayBlock)
  c1 = CachedArray(a*b)
  c2 = return_cache(unwrap_cached_array,c1)
  return c1, c2
end

function evaluate!(cache,::typeof(*),a::ArrayBlock,b::ArrayBlock)
  c1, c2 = cache
  setsize_op!(*,c1,a,b)
  c = evaluate!(c2,unwrap_cached_array,c1)
  mul!(c,a,b)
  c
end

# MulAddMap: Cached version of `mul!(d,a,b,α,β)`

"""
    struct MulAddMap{T} <: Map
      α::T
      β::T
    end

Map for the operation `a,b,c -> LinearAlgepra.mul!(copy(c),a,b,α,β)`
"""
struct MulAddMap{T} <: Map
  α::T
  β::T
end

function return_cache(k::MulAddMap,a,b,c)
  d = a*b+c
  CachedArray(d)
end

function evaluate!(cache,k::MulAddMap,a,b,c)
  setsize!(cache,size(c))
  d = cache.array
  copyto!(d,c)
  iszero(k.α) && isone(k.β) && return d
  mul!(d,a,b,k.α,k.β)
  d
end

function return_value(k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  x = return_value(*,a,b)
  return_value(+,x,c)
end

function return_cache(k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  c1 = CachedArray(a*b+c)
  c2 = return_cache(unwrap_cached_array,c1)
  (c1,c2)
end

function evaluate!(cache,k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  c1,c2 = cache
  setsize_op!(copy,c1,c)
  setsize_op!(*,c1,a,b)
  d = evaluate!(c2,unwrap_cached_array,c1)
  copyto!(d,c)
  iszero(k.α) && isone(k.β) && return d
  mul!(d,a,b,k.α,k.β)
  d
end

# Assembly Maps: AddEntriesMap and TouchEntriesMap

"""
    struct AddEntriesMap{F} <: Map

[`Map`](@ref) for [`add_entries!`](@ref) where `F` is the combine function.
"""
struct AddEntriesMap{F} <: Map
  combine::F
end

function evaluate!(cache,k::AddEntriesMap,A,v,i,j)
  add_entries!(k.combine,A,v,i,j)
end

function evaluate!(cache,k::AddEntriesMap,A,v,i)
  add_entries!(k.combine,A,v,i)
end

"""
    struct TouchEntriesMap <: Map

Dummy [`AddEntriesMap`](@ref) used to keep track of sparsity patterns in symbolic
assemblies, see the symbolic methods of [`SparseMatrixAssembler`](@ref Gridap.FESpaces.SparseMatrixAssembler).
"""
struct TouchEntriesMap <: Map end

function evaluate!(cache,k::TouchEntriesMap,A,v,i,j)
  add_entries!(+,A,nothing,i,j)
end

function evaluate!(cache,k::TouchEntriesMap,A,v,i)
  add_entries!(+,A,nothing,i)
end

for T in (:AddEntriesMap,:TouchEntriesMap)
  @eval begin

    function return_cache(k::$T,A,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
      qs = findall(v.touched)
      i, j = Tuple(first(qs))
      cij = return_cache(k,A,v.array[i,j],I.array[i],J.array[j])
      ni,nj = size(v.touched)
      cache = Matrix{typeof(cij)}(undef,ni,nj)
      for j in 1:nj
        for i in 1:ni
          if v.touched[i,j]
            cache[i,j] = return_cache(k,A,v.array[i,j],I.array[i],J.array[j])
          end
        end
      end
      cache
    end

    function evaluate!(cache, k::$T,A,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
      ni,nj = size(v.touched)
      for j in 1:nj
        for i in 1:ni
          if v.touched[i,j]
            evaluate!(cache[i,j],k,A,v.array[i,j],I.array[i],J.array[j])
          end
        end
      end
    end

    function return_cache(k::$T,A,v::VectorBlock,I::VectorBlock)
      qs = findall(v.touched)
      i = first(qs)
      ci = return_cache(k,A,v.array[i],I.array[i])
      ni = length(v.touched)
      cache = Vector{typeof(ci)}(undef,ni)
      for i in 1:ni
        if v.touched[i]
          cache[i] = return_cache(k,A,v.array[i],I.array[i])
        end
      end
      cache
    end

    function evaluate!(cache, k::$T,A,v::VectorBlock,I::VectorBlock)
      ni = length(v.touched)
      for i in 1:ni
        if v.touched[i]
          evaluate!(cache[i],k,A,v.array[i],I.array[i])
        end
      end
    end
  end

  for MT in (:MatrixBlock,:MatrixBlockView)
    Aij = (MT == :MatrixBlock) ? :(A.array[i,j]) : :(A[i,j])
    @eval begin

      function return_cache(k::$T,A::$MT,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
        qs = findall(v.touched)
        i, j = Tuple(first(qs))
        cij = return_cache(k,$Aij,v.array[i,j],I.array[i],J.array[j])
        ni,nj = size(v.touched)
        cache = Matrix{typeof(cij)}(undef,ni,nj)
        for j in 1:nj
          for i in 1:ni
            if v.touched[i,j]
              cache[i,j] = return_cache(k,$Aij,v.array[i,j],I.array[i],J.array[j])
            end
          end
        end
        cache
      end

      function evaluate!(cache,k::$T,A::$MT,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
        ni,nj = size(v.touched)
        for j in 1:nj
          for i in 1:ni
            if v.touched[i,j]
              evaluate!(cache[i,j],k,$Aij,v.array[i,j],I.array[i],J.array[j])
            end
          end
        end
      end

    end # @eval
  end # for MT

  for VT in (:VectorBlock,:VectorBlockView)
    Ai = (VT == :VectorBlock) ? :(A.array[i]) : :(A[i])
    @eval begin

      function return_cache(k::$T,A::$VT,v::VectorBlock,I::VectorBlock)
        qs = findall(v.touched)
        i = first(qs)
        ci = return_cache(k,$Ai,v.array[i],I.array[i])
        ni = length(v.touched)
        cache = Vector{typeof(ci)}(undef,ni)
        for i in 1:ni
          if v.touched[i]
            cache[i] = return_cache(k,$Ai,v.array[i],I.array[i])
          end
        end
        cache
      end

      function evaluate!(cache, k::$T,A::$VT,v::VectorBlock,I::VectorBlock)
        ni = length(v.touched)
        for i in 1:ni
          if v.touched[i]
            evaluate!(cache[i],k,$Ai,v.array[i],I.array[i])
          end
        end
      end

    end # @eval
  end # for VT
end # end for T

###########################################################################################
# Local solver maps

"""
    LocalSolveMap(; factorize! = lu!, pivot = RowMaximum())

A map for solving local linear systems `A x = B`, relying on a factorization method:
`x = ldiv!(factorize!(A,pivot),B)`.

- Given a matrix `A` and an array `B`, returns the array `x`.
- Given a block matrix `A` (test blocks x trial blocks) and a block array `B` (test blocks
  x input blocks, or a block vector for a single right-hand side), the blocks are flattened
  into a single system, which has to be square. Returns the solution as a block array
  (trial blocks x input blocks, or a block vector), which is reused between evaluations.

See the solution blocks returned by the other local maps in [`LocalStagedSolveMap`](@ref).
"""
struct LocalSolveMap{A,B} <: Map
  factorize! :: A
  pivot :: B
end

LocalSolveMap(; factorize! = lu!, pivot = RowMaximum()) = LocalSolveMap(factorize!, pivot)

Arrays.return_value(k::LocalSolveMap, matvec::Tuple) = return_value(k,matvec[1],matvec[2])
Arrays.return_cache(k::LocalSolveMap, matvec::Tuple) = return_cache(k,matvec[1],matvec[2])
Arrays.evaluate!(cache, k::LocalSolveMap, matvec::Tuple) = evaluate!(cache,k,matvec[1],matvec[2])

Arrays.return_value(k::LocalSolveMap, mat, vec) = similar(vec)

function Arrays.return_cache(k::LocalSolveMap, mat, vec)
  CachedArray(copy(mat)), CachedArray(copy(vec))
end

function Arrays.evaluate!(cache,k::LocalSolveMap, mat, vec)
  cmat, cvec = cache

  setsize!(cmat, size(mat))
  copy!(cmat.array, mat)
  f = k.factorize!(cmat.array,k.pivot;check=false)
  @check issuccess(f) "Factorization failed"

  setsize!(cvec, size(vec))
  ldiv!(cvec.array,f,vec)

  return cvec.array
end

# Block systems: a single stage with all the blocks, and all the trial blocks as output
# (see `LocalStagedSolveMap`)

function Arrays.return_value(k::LocalSolveMap, lhs::MatrixBlock, rhs::Union{MatrixBlock,VectorBlock})
  evaluate!(return_cache(k,lhs,rhs),k,lhs,rhs)
end

function Arrays.return_cache(k::LocalSolveMap, lhs::MatrixBlock, rhs::Union{MatrixBlock,VectorBlock})
  T = promote_type(eltype(eltype(lhs.array)), eltype(eltype(rhs.array)))
  offsets = _staged_offsets_cache(lhs, rhs)
  solved = fill(false, size(lhs, 2))
  X = CachedArray(zeros(T, 0, 0))
  stage = (CachedArray(zeros(T, 0, 0)), CachedArray(zeros(T, 0, 0)))
  out = _staged_output_cache(T, axes(lhs.array, 2), lhs, rhs)
  return offsets, solved, X, stage, out
end

function Arrays.evaluate!(cache, k::LocalSolveMap, lhs::MatrixBlock, rhs::Union{MatrixBlock,VectorBlock})
  @check size(lhs, 1) == size(rhs, 1)
  offsets, solved, cX, cstage, cout = cache
  ro, co, io = _staged_offsets!(offsets, lhs, rhs)

  setsize!(cX, (last(co), last(io)))
  X = cX.array
  fill!(solved, false)
  stage = (axes(lhs.array, 1), axes(lhs.array, 2))
  _staged_solve_stage!(k, stage, cstage, X, solved, offsets, lhs, rhs)

  return _staged_output!(cout, axes(lhs.array, 2), X, co, io)
end

###########################################################################################

"""
    LocalPenaltySolveMap(; factorize! = lu!, pivot = RowMaximum())

A map for solving local constrained linear systems, relying on a factorization method.

Given a left-hand-side 2x2 block matrix `mat` and a 2x2 block array (or a block vector of
length 2) of right-hand sides, returns the primal solution as a 1 x 2 block array (or a
block vector of length 1): the multiplier is not computed.

Each system is given by `A*[x_i; λ_i] = b_i`, where `A = [App, Aλp; Apλ, 0]` is the
augmented matrix, and `b_i = [Bp; Bλ]` is the right-hand side vector. The solution is
computed using a penalty method, as `x_i = ldiv!(factorize!(C,pivot),d_i)` with
`C = App + μT * Apλ * Aλp` and `d_i = Bp + μT * Apλ * Bλ`, where `μT` is a penalty parameter.
The penalty parameter μT is heuristically chosen as `μT = norm(App)/norm(Apλ*Aλp)`.

"""
struct LocalPenaltySolveMap{A,B} <: Map
  factorize! :: A
  pivot :: B
end

LocalPenaltySolveMap(; factorize! = lu!, pivot = RowMaximum()) = LocalPenaltySolveMap(factorize!, pivot)

Arrays.return_cache(k::LocalPenaltySolveMap, mats::Tuple) = return_cache(k, mats...)
Arrays.evaluate!(cache, k::LocalPenaltySolveMap, mats::Tuple) = evaluate!(cache, k, mats...)

# The solution blocks: the primal block, for each input block
Arrays.return_cache(k::LocalPenaltySolveMap, lhs, rhs) = _primal_solution_blocks(lhs, rhs)

function _primal_solution_blocks(lhs::MatrixBlock, rhs::VectorBlock)
  return ArrayBlock(Vector{eltype(rhs.array)}(undef, 1), fill(true, 1))
end

function _primal_solution_blocks(lhs::MatrixBlock, rhs::MatrixBlock)
  ni = size(rhs,2)
  return ArrayBlock(Matrix{eltype(rhs.array)}(undef, 1, ni), fill(true, 1, ni))
end

function Arrays.evaluate!(cache, k::LocalPenaltySolveMap, lhs, rhs::VectorBlock)
  @check size(lhs) == (2,2)
  @check size(rhs) == (2,)

  App, Aλp, Apλ, _ = get_array(lhs)
  Bp, Bλ = get_array(rhs)

  # μT = norm(App)/norm(Apλ*Aλp) is a heuristic choice for the penalty parameter
  if isone(size(Apλ,2))
    μT = tr(App)/norm(Apλ)^2 # Single constraint
  else
    μT = tr(App)/norm(Apλ*Aλp) # Multiple constraints
  end

  # App = App + μT * Apλ * Aλp
  mul!(App, Apλ, Aλp, μT, 1)

  # Bp = Bp + μT * Apλ * Bλ
  mul!(Bp, Apλ, Bλ, μT, 1)

  Ainv = k.factorize!(App,k.pivot;check=false)
  @check issuccess(Ainv) "Factorization failed"

  cache.array[1] = ldiv!(Ainv,Bp)
  return cache
end

function Arrays.evaluate!(cache, k::LocalPenaltySolveMap, lhs, rhs)
  @check size(lhs) == (2,2)
  @check size(rhs) == (2,2)

  App, Aλp, Apλ, _ = get_array(lhs)
  BpΩ, BλΩ, BpΓ, _ = get_array(rhs)

  # μT = norm(App)/norm(Apλ*Aλp) is a heuristic choice for the penalty parameter
  if isone(size(Apλ,2))
    μT = tr(App)/norm(Apλ)^2 # Single constraint
  else
    μT = tr(App)/norm(Apλ*Aλp) # Multiple constraints
  end

  # App = App + μT * Apλ * Aλp
  mul!(App, Apλ, Aλp, μT, 1)

  # BpΩ = BpΩ + μT * Apλ * BλΩ
  mul!(BpΩ, Apλ, BλΩ, μT, 1)

  Ainv = k.factorize!(App,k.pivot;check=false)
  @check issuccess(Ainv) "Factorization failed"

  cache.array[1,1] = ldiv!(Ainv,BpΩ)
  cache.array[1,2] = ldiv!(Ainv,BpΓ)
  return cache
end


"""
    LocalConstrainedSolveMap(; factorize! = lu!)

A map for solving local linear systems given as block arrays, with or without constraints.

Constrained case: given a 2x2 block matrix `lhs = [App Apλ; Aλp 0]` and a 2xN block array
`rhs`, whose columns are `[Bp_j; Bλ_j]`, computes the primal solutions `x_j` of

    [App Apλ; Aλp 0] * [x_j; λ_j] = [Bp_j; Bλ_j]

The multipliers `λ_j` are discarded. Blocks `Bλ_j` that are not touched are taken as zero.

Unconstrained case: given a 1x1 block matrix `lhs = [App]` and a 1xN block array `rhs`,
computes the solutions of `App * x_j = Bp_j`.

The primal solutions `x_j` are returned as a 1xN block array. If `rhs` is a block vector (a
single right-hand side), a block vector of length 1 is returned.

Unlike `LocalPenaltySolveMap`, the saddle-point system is solved exactly and
any number `N` of right-hand-side blocks is supported, each with its own constraint data.
"""
struct LocalConstrainedSolveMap{A} <: Map
  factorize! :: A
end

LocalConstrainedSolveMap(; factorize! = lu!) = LocalConstrainedSolveMap(factorize!)

Arrays.return_cache(k::LocalConstrainedSolveMap, mats::Tuple) = return_cache(k, mats...)
Arrays.evaluate!(cache, k::LocalConstrainedSolveMap, mats::Tuple) = evaluate!(cache, k, mats...)

Arrays.return_cache(k::LocalConstrainedSolveMap, lhs, rhs) = _primal_solution_blocks(lhs, rhs)

function Arrays.evaluate!(cache, k::LocalConstrainedSolveMap, lhs::MatrixBlock, rhs::VectorBlock)
  @check size(rhs,1) == size(lhs,1)
  F, np = _factorize_constrained_system(k, lhs)
  B = isone(size(rhs,1)) ? rhs.array[1] : vcat(rhs.array[1], rhs.array[2])
  cache.array[1] = (F \ B)[1:np]
  return cache
end

function Arrays.evaluate!(cache, k::LocalConstrainedSolveMap, lhs::MatrixBlock, rhs::MatrixBlock)
  @check size(rhs,1) == size(lhs,1)
  F, np = _factorize_constrained_system(k, lhs)
  nλ = size(F,1) - np
  constrained = (size(rhs,1) == 2)
  for j in axes(rhs.array,2)
    Bp = rhs.array[1,j]
    if constrained
      Bλ = rhs.touched[2,j] ? rhs.array[2,j] : zeros(eltype(Bp), nλ, size(Bp,2))
      B = vcat(Bp, Bλ)
    else
      B = Bp
    end
    x = F \ B
    cache.array[1,j] = x[1:np,:]
  end
  return cache
end

function _factorize_constrained_system(k::LocalConstrainedSolveMap, lhs::MatrixBlock)
  @check size(lhs) == (2,2) || size(lhs) == (1,1)
  App = lhs.array[1,1]
  np = size(App,1)
  if size(lhs) == (1,1)
    return k.factorize!(copy(App)), np
  end
  Aλp, Apλ = lhs.array[2,1], lhs.array[1,2]
  nλ = size(Aλp,1)
  K = [App Apλ; Aλp zeros(eltype(App), nλ, nλ)]
  return k.factorize!(K), np
end

"""
    LocalStagedSolveMap(stages; output = last trial block, factorize! = lu!, pivot = RowMaximum())

A map that solves local block systems `A x = B u` stage by stage, by block forward
substitution.

`A` is the block matrix assembled from the test fields (block rows) and the trial fields
(block columns) of a `LocalOperator`, and `B` holds one block column per input field.

`stages` is a tuple of `(rows, cols)` pairs, each a tuple of test block indices and a tuple
of trial block indices. The stages are solved in order: for stage `s`,

    x[cols_s] = A[rows_s, cols_s] \\ ( B[rows_s, :] - Σ_{t<s} A[rows_s, cols_t] * x[cols_t] )

where `A[rows_s, cols_s]` is flattened into a single matrix, which has to be square. The
blocks `A[rows_s, cols_t]` with `t > s` are ignored: the system is taken to be block lower
triangular in the stage ordering.

`output` is the trial block, or the tuple of trial blocks, that is returned. The other trial
blocks are intermediate unknowns: they are used by later stages and discarded.

The result is a block array with one row per block of `output`, in the order given, and one
column per input field (block column of `B`). If the right-hand side is a block vector (a
single function as input), it is a block vector with one entry per block of `output`. The
result is reused between evaluations.

All the local maps return the solution of a block system in this form. A `LocalOperator`
takes the rows, in order, as its output fields.

# Examples

- `LocalStagedSolveMap((((1,),(1,)),))`: a plain solve, as `LocalSolveMap`.
- `LocalStagedSolveMap((((1,2),(1,2)),); output = 1)`: a saddle-point system with a
  multiplier that is discarded.
- `LocalStagedSolveMap((((1,2),(1,)),))`: one unknown field tested against two families.
- `LocalStagedSolveMap((((1,),(1,)), ((2,3),(2,))))`: the composition of two operators. The
  first unknown is computed with the first test field. The second unknown is computed with
  the other two test fields, and its equations can involve the first unknown.
"""
struct LocalStagedSolveMap{S,O,A,B} <: Map
  stages :: S
  output :: O
  factorize! :: A
  pivot :: B
end

function LocalStagedSolveMap(
  stages; output = last(last(last(stages))), factorize! = lu!, pivot = RowMaximum()
)
  _stages = map(s -> (Tuple(s[1]), Tuple(s[2])), Tuple(stages))
  _output = Tuple(output)
  solved = reduce((a, b) -> (a..., b...), map(last, _stages))
  @check allunique(solved) "A trial block can only be solved in one stage"
  @check allunique(reduce((a, b) -> (a..., b...), map(first, _stages))) "A test block can only be used in one stage"
  @check all(in(solved), _output) "The output blocks have to be solved in some stage"
  LocalStagedSolveMap(_stages, _output, factorize!, pivot)
end

Arrays.return_cache(k::LocalStagedSolveMap, matvec::Tuple) = return_cache(k, matvec...)
Arrays.evaluate!(cache, k::LocalStagedSolveMap, matvec::Tuple) = evaluate!(cache, k, matvec...)

# Cache: block offsets of the test, trial and input blocks, which trial blocks are solved,
# the solution `X` of all the trial blocks (one column per input dof), the matrix and
# right-hand side of each stage, and the output blocks.
function Arrays.return_cache(k::LocalStagedSolveMap, lhs::MatrixBlock, rhs::Union{MatrixBlock,VectorBlock})
  T = promote_type(eltype(eltype(lhs.array)), eltype(eltype(rhs.array)))
  offsets = _staged_offsets_cache(lhs, rhs)
  solved = fill(false, size(lhs, 2))
  X = CachedArray(zeros(T, 0, 0))
  stages = map(_ -> (CachedArray(zeros(T, 0, 0)), CachedArray(zeros(T, 0, 0))), k.stages)
  out = _staged_output_cache(T, k.output, lhs, rhs)
  return offsets, solved, X, stages, out
end

# The solution blocks (output blocks x input blocks, or a vector of output blocks for a
# single right-hand side), and the arrays they point to.
function _staged_output_cache(::Type{T}, output, lhs::MatrixBlock, rhs::MatrixBlock) where T
  no, ni = length(output), size(rhs, 2)
  caches = [CachedArray(zeros(T, 0, 0)) for _ in 1:no, _ in 1:ni]
  return caches, ArrayBlock(Matrix{Matrix{T}}(undef, no, ni), fill(true, no, ni))
end
function _staged_output_cache(::Type{T}, output, lhs::MatrixBlock, rhs::VectorBlock) where T
  no = length(output)
  caches = [CachedArray(zeros(T, 0)) for _ in 1:no]
  return caches, ArrayBlock(Vector{Vector{T}}(undef, no), fill(true, no))
end

function Arrays.evaluate!(cache, k::LocalStagedSolveMap, lhs::MatrixBlock, rhs::Union{MatrixBlock,VectorBlock})
  @check size(lhs, 1) == size(rhs, 1)
  offsets, solved, cX, cstages, cout = cache
  ro, co, io = _staged_offsets!(offsets, lhs, rhs)

  setsize!(cX, (last(co), last(io)))
  X = cX.array
  fill!(solved, false)
  _staged_solve_stages!(k, k.stages, cstages, X, solved, offsets, lhs, rhs)

  return _staged_output!(cout, k.output, X, co, io)
end

# Copies the output blocks of the solution. For a single right-hand side, the blocks are
# vectors and there is a single input block.
function _staged_output!((caches, out), output, X, co, io)
  for b in axes(caches, 2), (o, j) in enumerate(output)
    c = caches[o, b]
    _staged_set_output_size!(c, co[j+1] - co[j], io[b+1] - io[b])
    copyto!(c.array, view(X, co[j]+1:co[j+1], io[b]+1:io[b+1]))
    out.array[o, b] = c.array
  end
  return out
end

# Recursion on the stages, which are tuples of different types
_staged_solve_stages!(k, ::Tuple{}, ::Tuple{}, args...) = nothing
function _staged_solve_stages!(k, stages::Tuple, cstages::Tuple, args...)
  _staged_solve_stage!(k, first(stages), first(cstages), args...)
  _staged_solve_stages!(k, Base.tail(stages), Base.tail(cstages), args...)
end

function _staged_solve_stage!(k, (rows, cols), (cF, cR), X, solved, (ro, co, io), lhs, rhs)
  n = _staged_num_entries(co, cols)
  @check _staged_num_entries(ro, rows) == n "Stage $((rows, cols)) is not square: $(_staged_num_entries(ro, rows)) equations, $n unknowns"
  setsize!(cF, (n, n))
  setsize!(cR, (n, last(io)))
  F, R = cF.array, cR.array
  fill!(F, zero(eltype(F)))
  fill!(R, zero(eltype(R)))

  oi = 0
  for i in rows
    I = oi+1:oi+ro[i+1]-ro[i]
    # Stage matrix
    oj = 0
    for j in cols
      J = oj+1:oj+co[j+1]-co[j]
      lhs.touched[i, j] && copyto!(view(F, I, J), lhs.array[i, j])
      oj = last(J)
    end
    # Right-hand side, minus the contributions of the unknowns solved in previous stages
    _staged_add_rhs!(R, I, io, rhs, i)
    for j in eachindex(solved)
      if solved[j] && lhs.touched[i, j]
        mul!(view(R, I, :), lhs.array[i, j], view(X, co[j]+1:co[j+1], :), -1, 1)
      end
    end
    oi = last(I)
  end

  f = k.factorize!(F, k.pivot; check = false)
  @check issuccess(f) "Factorization failed in stage $((rows, cols))"
  ldiv!(f, R)

  oj = 0
  for j in cols
    J = oj+1:oj+co[j+1]-co[j]
    copyto!(view(X, co[j]+1:co[j+1], :), view(R, J, :))
    solved[j] = true
    oj = last(J)
  end
  return nothing
end

_staged_set_output_size!(c::CachedVector, n, m) = setsize!(c, (n,))
_staged_set_output_size!(c::CachedMatrix, n, m) = setsize!(c, (n, m))
_staged_num_entries(offsets, blocks) = sum(b -> offsets[b+1] - offsets[b], blocks)

# Bases as input: one block column per input field
function _staged_add_rhs!(R, I, io, rhs::MatrixBlock, i)
  for b in axes(rhs.array, 2)
    rhs.touched[i, b] && copyto!(view(R, I, io[b]+1:io[b+1]), rhs.array[i, b])
  end
end

# A function as input: a single column
function _staged_add_rhs!(R, I, io, rhs::VectorBlock, i)
  rhs.touched[i] && copyto!(view(R, I, 1), rhs.array[i])
end

# Offsets of the test, trial and input blocks
function _staged_offsets_cache(lhs::MatrixBlock, rhs::Union{MatrixBlock,VectorBlock})
  (zeros(Int, size(lhs, 1) + 1), zeros(Int, size(lhs, 2) + 1), zeros(Int, size(rhs, 2) + 1))
end

function _staged_offsets!(offsets, lhs::MatrixBlock, rhs::Union{MatrixBlock,VectorBlock})
  ro, co, io = offsets
  _staged_block_offsets!(ro, lhs, 1)
  _staged_block_offsets!(co, lhs, 2)
  _staged_input_offsets!(io, rhs)
  return offsets
end

# Block offsets along dimension `dim`. The patch assembler returns every block allocated
# with its size, but cell-wise integration leaves the zero blocks untouched, so sizes are
# read from the touched ones.
function _staged_block_offsets!(offsets, m::MatrixBlock, dim)
  offsets[1] = 0
  for i in axes(m.array, dim)
    offsets[i+1] = offsets[i] + _staged_block_size(m, dim, i)
  end
  return offsets
end

function _staged_block_size(m::MatrixBlock, dim, i)
  for j in axes(m.array, 3 - dim)
    I = (dim == 1) ? CartesianIndex(i, j) : CartesianIndex(j, i)
    m.touched[I] && return size(m.array[I], dim)
  end
  @notimplemented "Cannot find the size of block $i (dimension $dim): it is empty"
end

_staged_input_offsets!(offsets, rhs::MatrixBlock) = _staged_block_offsets!(offsets, rhs, 2)
_staged_input_offsets!(offsets, rhs::VectorBlock) = (offsets[1] = 0; offsets[2] = 1; offsets)
