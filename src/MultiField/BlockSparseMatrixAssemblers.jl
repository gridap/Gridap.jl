
"""
    struct BlockSparseMatrixAssembler{R,C} <: FESpaces.SparseMatrixAssembler
        block_assemblers :: AbstractMatrix{<:Assembler}
    end

Block-wise sparse matrix assembler. This assembler is used to assemble 
block matrices, given an assembler for each block.

The block structure is given by the parameters `R` and `C`. Both are tuples containing 
the block structure of the rows and columns, respectively. These block structures are 
encoded in the following way: `R/C = (NB,SB,P)`

  - `NB` is the total number of blocks in the row/column direction.
  - `SB` is a tuple containing the number of fields in each block.
  - `P` is a tuple containing the permutation of the fields.

"""
struct BlockSparseMatrixAssembler{R,C,A} <: FESpaces.SparseMatrixAssembler
  block_assemblers :: A
  function BlockSparseMatrixAssembler{R,C}(
    block_assemblers::AbstractMatrix{<:Assembler}
  ) where {R,C}
    @assert isa(R,Tuple{<:Integer,<:Tuple,<:Tuple})
    @assert isa(C,Tuple{<:Integer,<:Tuple,<:Tuple})
    @assert R[1] == size(block_assemblers,1) == length(R[2])
    @assert C[1] == size(block_assemblers,2) == length(C[2])
    @assert sum(R[2]) == length(R[3])
    @assert sum(C[2]) == length(C[3])
    A = typeof(block_assemblers)
    return new{R,C,A}(block_assemblers)
  end
end

# Legacy constructor
function BlockSparseMatrixAssembler{NB,SB,P}(block_assemblers) where {NB,SB,P}
  R = (NB,SB,P)
  C = (NB,SB,P)
  return BlockSparseMatrixAssembler{R,C}(block_assemblers)
end

FESpaces.num_rows(a::BlockSparseMatrixAssembler) = sum(length,get_rows(a))
FESpaces.num_cols(a::BlockSparseMatrixAssembler) = sum(length,get_cols(a))

function FESpaces.get_rows(a::BlockSparseMatrixAssembler)
  assems = view(a.block_assemblers,:,1)
  return map(FESpaces.get_rows,assems)
end

function FESpaces.get_cols(a::BlockSparseMatrixAssembler)
  assems = view(a.block_assemblers,1,:)
  return map(FESpaces.get_cols,assems)
end

function FESpaces.get_assembly_strategy(a::BlockSparseMatrixAssembler)
  assems = a.block_assemblers
  strats = ArrayBlock(map(get_assembly_strategy,assems),fill(true,size(assems)))
  return expand_blocks(a,strats)
end

function FESpaces.get_matrix_builder(a::BlockSparseMatrixAssembler)
  assems = a.block_assemblers
  builders = ArrayBlock(map(get_matrix_builder,assems),fill(true,size(assems)))
  return expand_blocks(a,builders)
end

function FESpaces.get_vector_builder(a::BlockSparseMatrixAssembler)
  assems = view(a.block_assemblers,:,1)
  builders = ArrayBlock(map(get_vector_builder,assems),fill(true,length(assems)))
  return expand_blocks(a,builders)
end

function expand_blocks(a::BlockSparseMatrixAssembler,blocks::MatrixBlock)
  if has_trivial_blocks(a)
    return blocks
  end
  block_map = get_block_map(a)
  return ArrayBlockView(blocks,block_map)
end

function expand_blocks(a::BlockSparseMatrixAssembler,blocks::VectorBlock)
  if has_trivial_blocks(a)
    return blocks
  end
  block_map = map(idx -> CartesianIndex(idx[1]), view(get_block_map(a),:,1))
  return ArrayBlockView(blocks,block_map)
end

@inline function has_trivial_blocks(::BlockSparseMatrixAssembler{R,C}) where {R,C}
  return has_trivial_blocks(R...) && has_trivial_blocks(C...)
end

@inline get_block_map(::BlockSparseMatrixAssembler{R,C}) where {R,C} = get_block_map(R...,C...)

# Constructors

function BlockSparseMatrixAssembler(
  trial::MultiFieldFESpace,
  test::MultiFieldFESpace,
  matrix_builder,
  vector_builder,
  strategy=FESpaces.DefaultAssemblyStrategy()
)
  msg = "BlockSparseMatrixAssembler: trial and test spaces must have BlockMultiFieldStyle"
  @assert isa(MultiFieldStyle(trial),BlockMultiFieldStyle) msg
  @assert isa(MultiFieldStyle(test),BlockMultiFieldStyle) msg

  NBr, SBr, Pr = get_block_parameters(MultiFieldStyle(test))
  NBc, SBc, Pc = get_block_parameters(MultiFieldStyle(trial))

  # Count block rows/cols
  block_rows = [sum(num_free_dofs,test.spaces[r]) for r in get_block_ranges(NBr,SBr,Pr)]
  block_cols = [sum(num_free_dofs,trial.spaces[r]) for r in get_block_ranges(NBc,SBc,Pc)]

  # Create block assemblers
  block_idx = CartesianIndices((NBr,NBc))
  block_assemblers = map(block_idx) do idx
    rows = Base.OneTo(block_rows[idx[1]])
    cols = Base.OneTo(block_cols[idx[2]])
    FESpaces.GenericSparseMatrixAssembler(
      matrix_builder,vector_builder,rows,cols,strategy
    )
  end

  R, C = (NBr,SBr,Pr), (NBc,SBc,Pc)
  return BlockSparseMatrixAssembler{R,C}(block_assemblers)
end

function FESpaces.SparseMatrixAssembler(
  mat,vec,
  trial::MultiFieldFESpace{<:BlockMultiFieldStyle},
  test ::MultiFieldFESpace{<:BlockMultiFieldStyle},
  strategy::AssemblyStrategy=DefaultAssemblyStrategy()
)
  BlockSparseMatrixAssembler(
    trial,test,SparseMatrixBuilder(mat),ArrayBuilder(vec),strategy
  )
end

# BlockArrays extensions

function LinearAlgebra.fillstored!(a::BlockArray,v)
  foreach(ai->LinearAlgebra.fillstored!(ai,v),blocks(a))
  return a
end

# map cell ids

for T in [:MatrixBlock,:MatrixBlockView]
  @eval begin
    function FESpaces.map_cell_rows(strategy::$T{<:AssemblyStrategy},cell_ids,args...)
      strats = [strategy[i,1] for i in axes(strategy,1)]
      k = map(FESpaces.AssemblyStrategyMap{:rows},strats)
      return lazy_map(k,cell_ids,args...)
    end
    function FESpaces.map_cell_cols(strategy::$T{<:AssemblyStrategy},cell_ids,args...)
      strats = [strategy[1,i] for i in axes(strategy,2)]
      k = map(FESpaces.AssemblyStrategyMap{:cols},strats)
      return lazy_map(k,cell_ids,args...)
    end
    FESpaces.map_cell_rows(strategy::$T{FESpaces.DefaultAssemblyStrategy},cell_ids,args...) = cell_ids
    FESpaces.map_cell_cols(strategy::$T{FESpaces.DefaultAssemblyStrategy},cell_ids,args...) = cell_ids
  end
end

function Arrays.return_cache(k::Array{<:FESpaces.AssemblyStrategyMap},ids::ArrayBlock,args...)
  fi = testitem(ids)
  ki = testitem(k)
  ci = return_cache(ki,fi,args...)
  gi = evaluate!(ci,ki,fi,args...)
  b = Array{typeof(ci),ndims(ids)}(undef,size(ids))
  for i in eachindex(ids.array)
    if ids.touched[i]
      ki = return_cache(k[i],ids.array[i],args...)
      b[i] = return_cache(k[i],ids.array[i],args...)
    end
  end
  array = Array{typeof(gi),ndims(ids)}(undef,size(ids))
  ArrayBlock(array,ids.touched), b
end

function Arrays.evaluate!(cache,k::Array{<:FESpaces.AssemblyStrategyMap},ids::ArrayBlock,args...)
  a,b = cache
  for i in eachindex(ids.array)
    if ids.touched[i]
      a.array[i] = evaluate!(b[i],k[i],ids.array[i],args...)
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

function FESpaces.assemble_matrix_and_vector_add!(
  A::AbstractBlockMatrix, b::AbstractBlockVector, a::BlockSparseMatrixAssembler, data
)
  m1 = ArrayBlock(blocks(A),fill(true,blocksize(A)))
  m2 = expand_blocks(a,m1)
  b1 = ArrayBlock(blocks(b),fill(true,blocksize(b)))
  b2 = expand_blocks(a,b1)
  FESpaces.assemble_matrix_and_vector_add!(m2,b2,a,data)
end
