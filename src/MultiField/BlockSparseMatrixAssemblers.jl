
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

function FESpaces.get_assembly_strategy(a::BlockSparseMatrixAssembler)
  assems = a.block_assemblers
  return ArrayBlock(map(get_assembly_strategy,assems),fill(true,size(assems)))
end

function FESpaces.get_matrix_builder(a::BlockSparseMatrixAssembler)
  assems = a.block_assemblers
  return ArrayBlock(map(get_matrix_builder,assems),fill(true,size(assems)))
end

function FESpaces.get_vector_builder(a::BlockSparseMatrixAssembler)
  assems = diag(a.block_assemblers)
  return ArrayBlock(map(get_vector_builder,assems),fill(true,length(assems)))
end

function get_block_map(::BlockSparseMatrixAssembler{NB,NV,SB,P}) where {NB,NV,SB,P}
  #! TODO: Include permutation P into the block map
  ptrs = [1,SB...]
  length_to_ptrs!(ptrs)

  block_map = Matrix{CartesianIndex{2}}(undef,NV,NV)
  for I in CartesianIndices((NB,NB))
    i_range = ptrs[I[1]]:ptrs[I[1]+1]-1
    j_range = ptrs[I[2]]:ptrs[I[2]+1]-1
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

function BlockSparseMatrixAssembler(trial::MultiFieldFESpace{<:BlockMultiFieldStyle{NB,SB,P}},
                                    test::MultiFieldFESpace{<:BlockMultiFieldStyle{NB,SB,P}},
                                    matrix_builder,
                                    vector_builder,
                                    strategy=FESpaces.DefaultAssemblyStrategy()) where {NB,SB,P}

  block_idx = CartesianIndices((length(test),length(trial)))
  block_assemblers = map(block_idx) do idx
    block_rows = get_free_dof_ids(test[idx[1]])
    block_cols = get_free_dof_ids(trial[idx[2]])
    FESpaces.GenericSparseMatrixAssembler(matrix_builder,vector_builder,block_rows,block_cols,strategy)
  end

  NV = length(test.spaces)
  return BlockSparseMatrixAssembler{NB,NV,SB,P}(block_assemblers)
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

struct ArrayBlockView{A,N,M}
  array::ArrayBlock{A,M}
  block_map::Array{CartesianIndex{M},N}
end

Base.view(a::ArrayBlock{A,M},b::Array{CartesianIndex{M},N}) where {A,M,N} = ArrayBlockView(a,b)
const MatrixBlockView{A} = ArrayBlockView{A,2,2} where A
const VectorBlockView{A} = ArrayBlockView{A,1,1} where A

Base.size(a::ArrayBlockView) = size(a.block_map)
Base.length(b::ArrayBlockView) = length(b.block_map)
Base.eltype(::Type{<:ArrayBlockView{A}}) where A = A
Base.eltype(::ArrayBlockView{A}) where A = A
Base.ndims(::ArrayBlockView{A,N}) where {A,N} = N
Base.ndims(::Type{ArrayBlockView{A,N}}) where {A,N} = N
Base.getindex(b::ArrayBlockView,i...) = getindex(b.array,b.block_map[i...])
Base.setindex!(b::ArrayBlockView,v,i...) = setindex!(b.array,v,b.block_map[i...])

Base.copy(a::ArrayBlockView) = ArrayBlockView(copy(a.array),copy(a.block_map))
Base.eachindex(a::ArrayBlockView) = eachindex(a.block_map)

function Base.show(io::IO,o::ArrayBlockView)
  print(io,"ArrayBlockView($(o.array), $(o.block_map))")
end

LinearAlgebra.diag(a::MatrixBlock) = view(a.array,diag(CartesianIndices(a.array)))
LinearAlgebra.diag(a::MatrixBlockView) = view(a.array.array, diag(a.block_map))

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

  for MT in (:MatrixBlock,:MatrixBlockView,:BlockMatrix)
    Aij = (MT == :BlockMatrix) ? :(A.blocks[i,j]) : :(A[i,j])
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

  for VT in (:VectorBlock,:VectorBlockView,:BlockVector)
    Ai = (VT == :BlockVector) ? :(A.blocks[i]) : :(A[i])
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
