

struct BlockSparseMatrixAssembler{A} <: FESpaces.SparseMatrixAssembler
  block_assemblers :: AbstractMatrix{A}
end

FESpaces.num_rows(a::BlockSparseMatrixAssembler) = sum(map(length,get_rows(a)))
FESpaces.num_cols(a::BlockSparseMatrixAssembler) = sum(map(length,get_cols(a)))

function FESpaces.get_rows(a::BlockSparseMatrixAssembler)
  row_assemblers = a.block_assemblers[:,1]
  return map(FESpaces.get_rows,row_assemblers)
end

function FESpaces.get_cols(a::BlockSparseMatrixAssembler)
  col_assemblers = a.block_assemblers[1,:]
  return map(FESpaces.get_cols,col_assemblers)
end

function FESpaces.get_assembly_strategy(a::BlockSparseMatrixAssembler)
  return get_assembly_strategy(first(a.block_assemblers))
end

function FESpaces.get_matrix_builder(a::BlockSparseMatrixAssembler)
  assems = a.block_assemblers
  return ArrayBlock(map(get_matrix_builder,assems),fill(true,size(assems)))
end

function FESpaces.get_vector_builder(a::BlockSparseMatrixAssembler)
  assems = diag(a.block_assemblers)
  return ArrayBlock(map(get_vector_builder,assems),fill(true,length(assems)))
end

# Constructors

function BlockSparseMatrixAssembler(trial::MultiFieldFESpace{<:MS},
                                    test::MultiFieldFESpace{<:MS},
                                    matrix_builder,
                                    vector_builder,
                                    strategy=FESpaces.DefaultAssemblyStrategy()) where MS
  msg = "Block assembly is only allowed for BlockMultiFieldStyle."
  @check (MS <: BlockMultiFieldStyle) msg

  block_idx = CartesianIndices((length(test),length(trial)))
  block_assemblers = map(block_idx) do idx
    block_rows = get_free_dof_ids(test[idx[1]])
    block_cols = get_free_dof_ids(trial[idx[2]])
    FESpaces.GenericSparseMatrixAssembler(matrix_builder,vector_builder,block_rows,block_cols,strategy)
  end

  return BlockSparseMatrixAssembler(block_assemblers)
end

function FESpaces.SparseMatrixAssembler(mat,
                                        vec,
                                        trial::MultiFieldFESpace{<:BlockMultiFieldStyle},
                                        test::MultiFieldFESpace{<:BlockMultiFieldStyle},
                                        strategy::AssemblyStrategy=DefaultAssemblyStrategy())
  return BlockSparseMatrixAssembler(trial,test,SparseMatrixBuilder(mat),ArrayBuilder(vec),strategy)
end

# BlockArrays extensions

function LinearAlgebra.fillstored!(a::BlockArray,v)
  map(ai->LinearAlgebra.fillstored!(ai,v),blocks(a))
end

# nnz counters and allocators

Algebra.LoopStyle(a::ArrayBlock) = Algebra.LoopStyle(first(a.array))

function Algebra.nz_counter(builder::MatrixBlock,axs)
  s = size(builder.array)
  rows = axs[1]
  cols = axs[2]
  counters = [nz_counter(builder.array[i,j],(rows[i],cols[j])) for i in 1:s[1], j in 1:s[2]]
  return ArrayBlock(counters,fill(true,size(counters)))
end

function Algebra.nz_counter(builder::VectorBlock,axs)
  s = size(builder.array)
  rows = axs[1]
  counters = [nz_counter(builder.array[i],(rows[i],)) for i in 1:s[1]]
  return ArrayBlock(counters,fill(true,size(counters)))
end

function Algebra.nz_allocation(a::ArrayBlock)
  array = map(Algebra.nz_allocation,a.array)
  return ArrayBlock(array,a.touched)
end

function Algebra.create_from_nz(a::ArrayBlock)
  array = map(Algebra.create_from_nz,a.array)
  return mortar(array)
end

# AddEntriesMap and TouchEntriesMap

for T in (:AddEntriesMap,:TouchEntriesMap)

  for MT in (:MatrixBlock,BlockMatrix)
    Aij = (MT == :MatrixBlock) ? :(A.array[i,j]) : :(A.blocks[i,j])
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

      function Fields.evaluate!(cache, k::$T,A::$MT,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
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

  for VT in (:VectorBlock,BlockVector)
    Ai = (VT == :VectorBlock) ? :(A.array[i]) : :(A.blocks[i])
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
