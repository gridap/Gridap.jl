
struct BlockSparseMatrixAssembler <: FESpaces.SparseMatrixAssembler
  glob_assembler   :: SparseMatrixAssembler
  block_assemblers :: AbstractArray{<:SparseMatrixAssembler}

  function BlockSparseMatrixAssembler(trial::MultiFieldFESpace{<:MS},
                                      test::MultiFieldFESpace{<:MS},
                                      matrix_builder,
                                      vector_builder,
                                      strategy=FESpaces.DefaultAssemblyStrategy()) where MS
    msg = "Block assembly is only allowed for BlockMultiFieldStyle."
    @check (MS <: BlockMultiFieldStyle) msg

    # Regular global assembler
    rows = get_free_dof_ids(test)
    cols = get_free_dof_ids(trial)
    glob_assembler = FESpaces.GenericSparseMatrixAssembler(matrix_builder,
                                                           vector_builder,
                                                           rows,
                                                           cols,
                                                           strategy)

    # Block assemblers
    nblocks = (length(test),length(trial))
    block_assemblers = Matrix{SparseMatrixAssembler}(undef,nblocks)
    for i in 1:nblocks[1], j in 1:nblocks[2]
      block_rows = get_free_dof_ids(test[i])
      block_cols = get_free_dof_ids(trial[j])
      block_assemblers[i,j] = FESpaces.GenericSparseMatrixAssembler(matrix_builder,
                                                                    vector_builder,
                                                                    block_rows,
                                                                    block_cols,
                                                                    strategy)
    end
    
    return new{}(glob_assembler,block_assemblers)
  end
end

function FESpaces.SparseMatrixAssembler(
  mat,
  vec,
  trial::MultiFieldFESpace{<:BlockMultiFieldStyle},
  test::MultiFieldFESpace{<:BlockMultiFieldStyle},
  strategy::AssemblyStrategy=DefaultAssemblyStrategy())
  return BlockSparseMatrixAssembler(trial,test,SparseMatrixBuilder(mat),ArrayBuilder(vec),strategy)
end

for fun in [:get_rows,:get_cols,:get_matrix_builder,:get_vector_builder,:get_assembly_strategy]
  @eval begin
    function FESpaces.$fun(a::BlockSparseMatrixAssembler)
      $fun(a.glob_assembler)
    end
  end
end

function allocate_block_vector(a::BlockSparseMatrixAssembler)
  vec_type = get_vector_builder(a.glob_assembler).array_type
  rows = get_rows(a.glob_assembler)
  r = rows.lasts .- [0,rows.lasts[1:end-1]...]
  BlockArray(undef_blocks,vec_type,r)
end

function allocate_block_matrix(a::BlockSparseMatrixAssembler)
  mat_type = get_matrix_builder(a.glob_assembler).matrix_type
  rows = get_rows(a.glob_assembler)
  cols = get_cols(a.glob_assembler)
  r = rows.lasts .- [0,rows.lasts[1:end-1]...]
  c = cols.lasts .- [0,cols.lasts[1:end-1]...]
  BlockArray(undef_blocks,mat_type,r,c)
end

function select_block_matdata(matdata,i::Integer,j::Integer)
  (map(data->_select_block_data(data,i,j),matdata[1]),
   map(data->_select_block_data(data,i),matdata[2]),
   map(data->_select_block_data(data,j),matdata[3]))
end

function select_block_vecdata(vecdata,j::Integer)
  (map(data->_select_block_data(data,j),vecdata[1]),
   map(data->_select_block_data(data,j),vecdata[2]))
end

function _select_block_data(data,idx::Integer...)
  n = length(data)
  _idx = map(i->Fill(i,n),idx)
  return lazy_map(getindex,data,_idx...)
end

"""
  TODO: We need to detect inactive blocks and avoid assembling them.
  Otherwise, we allocate unnecessary memory.
"""

# Vector assembly 

function FESpaces.assemble_vector(a::BlockSparseMatrixAssembler,vecdata)
  v = allocate_block_vector(a)
  block_assemblers = a.block_assemblers
  for j in 1:blocksize(v,1)
    _a = block_assemblers[j,1]
    _vecdata = select_block_vecdata(vecdata,j)
    v[Block(j)] = assemble_vector(_a,_vecdata)
  end
  return v
end

function FESpaces.allocate_vector(a::BlockSparseMatrixAssembler,vecdata)
  v = allocate_block_vector(a)
  block_assemblers = a.block_assemblers
  for j in 1:blocksize(v,1)
    _a = block_assemblers[j,1]
    _vecdata = select_block_vecdata(vecdata,j)
    v[Block(i)] = allocate_vector(_a,_vecdata)
  end
  return v
end

function FESpaces.assemble_vector_add!(b::BlockVector,a::BlockSparseMatrixAssembler,vecdata)
  block_assemblers = a.block_assemblers
  for j in 1:blocksize(v,1)
    _a = block_assemblers[j,1]
    _vecdata = select_block_vecdata(vecdata,j)
    bj = b[Block(j)]
    assemble_vector_add!(bj,_a,_vecdata)
  end
end

# Matrix assembly 

function FESpaces.assemble_matrix(a::BlockSparseMatrixAssembler,matdata)
  m = allocate_block_matrix(a)
  block_assemblers = a.block_assemblers
  for i in 1:blocksize(m,1)
    for j in 1:blocksize(m,2)
      _a = block_assemblers[i,j]
      _matdata = select_block_matdata(matdata,i,j)
      m[Block(i,j)] = assemble_matrix(_a,_matdata)
    end
  end
  return m
end

function FESpaces.allocate_matrix(a::BlockSparseMatrixAssembler,matdata)
  m = allocate_block_matrix(a)
  block_assemblers = a.block_assemblers
  for i in 1:blocksize(m,1)
    for j in 1:blocksize(m,2)
      _a = block_assemblers[i,j]
      _matdata = select_block_matdata(matdata,i,j)
      m[Block(i,j)] = allocate_matrix(_a,_matdata)
    end
  end
  return m
end

function FESpaces.assemble_matrix!(mat::BlockMatrix,a::BlockSparseMatrixAssembler,matdata)
  for (i,j) in blockaxes(mat)
    LinearAlgebra.fillstored!(mat[i,j],zero(eltype(mat)))
  end
  assemble_matrix_add!(mat,a,matdata)
end

function FESpaces.assemble_matrix_add!(mat::BlockMatrix,a::BlockSparseMatrixAssembler,matdata)
  block_assemblers = a.block_assemblers
  for i in 1:blocksize(m,1)
    for j in 1:blocksize(m,2)
      _a = block_assemblers[i,j]
      _matdata = select_block_matdata(matdata,i,j)
      mat_ij = mat[Block(i,j)]
      assemble_matrix_add!(mat_ij,_a,_matdata)
    end
  end
end

# Matrix and vector Assembly

"""
function FESpaces.allocate_matrix_and_vector(a::BlockSparseMatrixAssembler,data)
  
end

function FESpaces.assemble_matrix_and_vector!(A,b,a::BlockSparseMatrixAssembler, data)

end

function FESpaces.assemble_matrix_and_vector_add!(A,b,a::BlockSparseMatrixAssembler,data)

end

function FESpaces.assemble_matrix_and_vector(a::BlockSparseMatrixAssembler,data)

end
"""

