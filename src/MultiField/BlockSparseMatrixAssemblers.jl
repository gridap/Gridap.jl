
struct BlockSparseMatrixAssembler{NB,NV,SB,P,A} <: FESpaces.SparseMatrixAssembler
  block_assemblers :: AbstractMatrix{A}
  function BlockSparseMatrixAssembler{NB,NV,SB,P}(block_assemblers) where {NB,NV,SB,P}
    A = eltype(block_assemblers)
    return new{NB,NV,SB,P,A}(block_assemblers)
  end
end

FESpaces.num_rows(a::BlockSparseMatrixAssembler) = sum(map(length,get_rows(a)))
FESpaces.num_cols(a::BlockSparseMatrixAssembler) = sum(map(length,get_cols(a)))

function FESpaces.get_rows(a::BlockSparseMatrixAssembler)
  return map(FESpaces.get_rows,diag(a.block_assemblers))
end

function FESpaces.get_cols(a::BlockSparseMatrixAssembler)
  return map(FESpaces.get_cols,diag(a.block_assemblers))
end

function FESpaces.get_assembly_strategy(a::BlockSparseMatrixAssembler{NB,NV}) where {NB,NV}
  assems = a.block_assemblers
  strats = ArrayBlock(map(get_assembly_strategy,assems),fill(true,size(assems)))
  return expand_blocks(a,strats)
end

function FESpaces.get_matrix_builder(a::BlockSparseMatrixAssembler{NB,NV}) where {NB,NV}
  assems = a.block_assemblers
  builders = ArrayBlock(map(get_matrix_builder,assems),fill(true,size(assems)))
  return expand_blocks(a,builders)
end

function FESpaces.get_vector_builder(a::BlockSparseMatrixAssembler{NB,NV}) where {NB,NV}
  assems = diag(a.block_assemblers)
  builders = ArrayBlock(map(get_vector_builder,assems),fill(true,length(assems)))
  return expand_blocks(a,builders)
end

function expand_blocks(a::BlockSparseMatrixAssembler{NB,NV,SB,P},blocks::MatrixBlock) where {NB,NV,SB,P}
  if (NB == NV) && all(x -> x[1] == x[2], enumerate(P))
    return blocks
  end
  block_map = get_block_map(a)
  return ArrayBlockView(blocks,block_map)
end

function expand_blocks(a::BlockSparseMatrixAssembler{NB,NV,SB,P},blocks::VectorBlock) where {NB,NV,SB,P}
  if (NB == NV) && all(x -> x[1] == x[2], enumerate(P))
    return blocks
  end
  block_map = map(idx -> CartesianIndex(idx[1]), diag(get_block_map(a)))
  return ArrayBlockView(blocks,block_map)
end

function get_block_ranges(NB::Integer,SB,P)
  ptrs = [1,SB...]
  length_to_ptrs!(ptrs)
  var_perm = [P...]
  return map(i-> var_perm[ptrs[i]:ptrs[i+1]-1], 1:NB)
end

function get_block_map(::BlockSparseMatrixAssembler{NB,NV,SB,P}) where {NB,NV,SB,P}
  ranges = get_block_ranges(NB,SB,P)
  block_map = Matrix{CartesianIndex{2}}(undef,NV,NV)
  for I in CartesianIndices((NB,NB))
    i_range = ranges[I[1]]
    j_range = ranges[I[2]]
    for i in i_range, j in j_range
      block_map[i,j] = I
    end
  end
  return block_map
end

# Constructors

function BlockSparseMatrixAssembler(trial::MultiFieldFESpace,
                                    test::MultiFieldFESpace,
                                    matrix_builder,
                                    vector_builder,
                                    strategy=FESpaces.DefaultAssemblyStrategy())
  msg = "Block assembly is only allowed for BlockMultiFieldStyle."
  @notimplemented msg
end

function BlockSparseMatrixAssembler(::BlockMultiFieldStyle{NB,SB,P},
                                    trial,
                                    test,
                                    matrix_builder,
                                    vector_builder,
                                    strategy=FESpaces.DefaultAssemblyStrategy()) where {NB,SB,P}

  # Count block rows/cols
  NV = length(test.spaces)
  block_ranges = get_block_ranges(NB,SB,P)
  block_rows = map(range->sum(map(num_free_dofs,test.spaces[range])),block_ranges)
  block_cols = map(range->sum(map(num_free_dofs,trial.spaces[range])),block_ranges)

  # Create block assemblers
  block_idx = CartesianIndices((NB,NB))
  block_assemblers = map(block_idx) do idx
    rows = Base.OneTo(block_rows[idx[1]])
    cols = Base.OneTo(block_cols[idx[2]])
    FESpaces.GenericSparseMatrixAssembler(matrix_builder,vector_builder,rows,cols,strategy)
  end

  return BlockSparseMatrixAssembler{NB,NV,SB,P}(block_assemblers)
end

function FESpaces.SparseMatrixAssembler(mat,vec,
                                        trial::MultiFieldFESpace{MS},
                                        test ::MultiFieldFESpace{MS},
                                        strategy::AssemblyStrategy=DefaultAssemblyStrategy()
                                       ) where MS <: BlockMultiFieldStyle
  mfs = MultiFieldStyle(test)
  return BlockSparseMatrixAssembler(mfs,trial,test,SparseMatrixBuilder(mat),ArrayBuilder(vec),strategy)
end

# BlockArrays extensions

function LinearAlgebra.fillstored!(a::BlockArray,v)
  map(ai->LinearAlgebra.fillstored!(ai,v),blocks(a))
end

# map cell ids

for T in [:MatrixBlock,:MatrixBlockView]
  @eval begin
    function FESpaces.map_cell_rows(strategy::$T{<:AssemblyStrategy},cell_ids)
      strats = diag(strategy)
      k = map(FESpaces.AssemblyStrategyMap{:rows},strats)
      return lazy_map(k,cell_ids)
    end
    function FESpaces.map_cell_cols(strategy::$T{<:AssemblyStrategy},cell_ids)
      strats = diag(strategy)
      k = map(FESpaces.AssemblyStrategyMap{:cols},strats)
      return lazy_map(k,cell_ids)
    end
  end
end

function Arrays.return_cache(k::Array{<:FESpaces.AssemblyStrategyMap},ids::ArrayBlock)
  fi = testitem(ids)
  ki = testitem(k)
  ci = return_cache(ki,fi)
  gi = evaluate!(ci,ki,fi)
  b = Array{typeof(ci),ndims(ids)}(undef,size(ids))
  for i in eachindex(ids.array)
    if ids.touched[i]
      ki = return_cache(k[i],ids.array[i])
      b[i] = return_cache(k[i],ids.array[i])
    end
  end
  array = Array{typeof(gi),ndims(ids)}(undef,size(ids))
  ArrayBlock(array,ids.touched), b
end

function Arrays.evaluate!(cache,k::Array{<:FESpaces.AssemblyStrategyMap},ids::ArrayBlock)
  a,b = cache
  for i in eachindex(ids.array)
    if ids.touched[i]
      a.array[i] = evaluate!(b[i],k[i],ids.array[i])
    end
  end
  a
end

# nnz counters and allocators

Algebra.LoopStyle(a::ArrayBlock) = Algebra.LoopStyle(first(a.array))
Algebra.LoopStyle(a::ArrayBlockView) = Algebra.LoopStyle(a.array)

function Algebra.nz_counter(builder::MatrixBlock,axs)
  s = size(builder)
  rows = axs[1]
  cols = axs[2]
  counters = [nz_counter(builder.array[i,j],(rows[i],cols[j])) for i in 1:s[1], j in 1:s[2]]
  return ArrayBlock(counters,fill(true,size(counters)))
end

function Algebra.nz_counter(builder::VectorBlock,axs)
  s = size(builder)
  rows = axs[1]
  counters = [nz_counter(builder.array[i],(rows[i],)) for i in 1:s[1]]
  return ArrayBlock(counters,fill(true,size(counters)))
end

function Algebra.nz_counter(builder::ArrayBlockView,axs)
  ArrayBlockView(nz_counter(builder.array,axs),builder.block_map)
end

function Algebra.nz_allocation(a::ArrayBlock)
  array = map(Algebra.nz_allocation,a.array)
  return ArrayBlock(array,a.touched)
end

function Algebra.nz_allocation(a::ArrayBlockView)
  ArrayBlockView(nz_allocation(a.array),a.block_map)
end

function Algebra.create_from_nz(a::ArrayBlock)
  array = map(Algebra.create_from_nz,a.array)
  return mortar(array)
end

Algebra.create_from_nz(a::ArrayBlockView) = create_from_nz(a.array)

# AddEntriesMap and TouchEntriesMap

for T in (:AddEntriesMap,:TouchEntriesMap)

  for MT in (:MatrixBlock,:MatrixBlockView)
    Aij = (MT == :MatrixBlock) ? :(A.array[i,j]) : :(A[i,j])
    @eval begin

      function Fields.return_cache(k::$T,A::$MT,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
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

      function Fields.evaluate!(cache,k::$T,A::$MT,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
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

      function Fields.return_cache(k::$T,A::$VT,v::VectorBlock,I::VectorBlock)
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

      function Fields.evaluate!(cache, k::$T,A::$VT,v::VectorBlock,I::VectorBlock)
        ni = length(v.touched)
        for i in 1:ni
          if v.touched[i]
            evaluate!(cache[i],k,$Ai,v.array[i],I.array[i])
          end
        end
      end

    end # @eval
  end # for VT
end

# In place assembly modifications (for dispatching)
# We convert from BlockArray to ArrayBlock to be able to expland the blocks. 
# After assembly we convert back to BlockArray automatically.
function FESpaces.assemble_vector_add!(b::AbstractBlockVector,a::BlockSparseMatrixAssembler,vecdata)
  b1 = ArrayBlock(blocks(b),fill(true,blocksize(b)))
  b2 = expand_blocks(a,b1)
  FESpaces.assemble_vector_add!(b2,a,vecdata)
end

function FESpaces.assemble_matrix_add!(mat::AbstractBlockMatrix,a::BlockSparseMatrixAssembler,matdata)
  m1 = ArrayBlock(blocks(mat),fill(true,blocksize(mat)))
  m2 = expand_blocks(a,m1)
  FESpaces.assemble_matrix_add!(m2,a,matdata)
end

function FESpaces.assemble_matrix_and_vector_add!(A::AbstractBlockMatrix,
                                                  b::AbstractBlockVector,
                                                  a::BlockSparseMatrixAssembler,
                                                  data)
  m1 = ArrayBlock(blocks(A),fill(true,blocksize(A)))
  m2 = expand_blocks(a,m1)
  b1 = ArrayBlock(blocks(b),fill(true,blocksize(b)))
  b2 = expand_blocks(a,b1)
  FESpaces.assemble_matrix_and_vector_add!(m2,b2,a,data)
end