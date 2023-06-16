
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

function select_block_matvecdata(matvecdata,i::Integer,j::Integer)
  (map(data->lazy_map(_select_block_matvecdata{i,j}(),data),matvecdata[1]),
   map(data->_select_block_data(data,i),matvecdata[2]),
   map(data->_select_block_data(data,j),matvecdata[3]))
end

function _select_block_data(data,idx::Integer...)
  n = length(data)
  _idx = map(i->Fill(i,n),idx)
  return lazy_map(getindex,data,_idx...)
end

struct _select_block_matvecdata{I,J} <: Map end

function Arrays.evaluate!(cache,::_select_block_matvecdata{I,J},data) where {I,J}
  matdata, vecdata = data
  return matdata[I,J]
end

function Arrays.evaluate!(cache,::_select_block_matvecdata{I,I},data) where {I}
  matdata, vecdata = data
  return matdata[I,I], vecdata[I]
end

function combine_matdata(data1,data2)
  w1,r1,c1 = data1
  w2,r2,c2 = data2
  return (vcat(w1,w2),vcat(r1,r2),vcat(c1,c2))
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
  for i in 1:blocksize(v,1)
    _a = block_assemblers[i,1]
    _vecdata = select_block_vecdata(vecdata,i)
    v[Block(i)] = allocate_vector(_a,_vecdata)
  end
  return v
end

function FESpaces.assemble_vector_add!(b::BlockVector,a::BlockSparseMatrixAssembler,vecdata)
  block_assemblers = a.block_assemblers
  for j in 1:blocksize(b,1)
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
  for i in 1:blocksize(mat,1)
    for j in 1:blocksize(mat,2)
      _a = block_assemblers[i,j]
      _matdata = select_block_matdata(matdata,i,j)
      mat_ij = mat[Block(i,j)]
      assemble_matrix_add!(mat_ij,_a,_matdata)
    end
  end
end

# Matrix and vector Assembly

function FESpaces.allocate_matrix_and_vector(a::BlockSparseMatrixAssembler,data)
  matvecdata, matdata, vecdata = data
  m = allocate_block_matrix(a)
  v = allocate_block_vector(a)
  block_assemblers = a.block_assemblers
  for i in 1:blocksize(m,1)
    for j in 1:blocksize(m,2)
      _a = block_assemblers[i,j]
      _matvecdata = select_block_matvecdata(matvecdata,i,j)
      _matdata = select_block_matdata(matdata,i,j)
      if i == j # Diagonal blocks
        _vecdata = select_block_vecdata(vecdata,j)
        m[Block(i,j)], v[Block(i)] = allocate_matrix_and_vector(_a,(_matvecdata,_matdata,_vecdata))
      else      # Off-diagonal blocks
        __matdata = combine_matdata(_matvecdata,_matdata) 
        m[Block(i,j)] = allocate_matrix(_a,__matdata)
      end
    end
  end
  return m, v
end

function FESpaces.assemble_matrix_and_vector!(A::BlockMatrix,b::BlockVector,a::BlockSparseMatrixAssembler,data)
  for (i,j) in blockaxes(A)
    LinearAlgebra.fillstored!(A[i,j],zero(eltype(A)))
  end
  fill!(b,zero(eltype(b)))
  assemble_matrix_and_vector_add!(A,b,a,data)
end

function FESpaces.assemble_matrix_and_vector_add!(A::BlockMatrix,b::BlockVector,a::BlockSparseMatrixAssembler,data)
  matvecdata, matdata, vecdata = data
  block_assemblers = a.block_assemblers
  for i in 1:blocksize(A,1)
    for j in 1:blocksize(A,2)
      _a = block_assemblers[i,j]
      _matvecdata = select_block_matvecdata(matvecdata,i,j)
      _matdata = select_block_matdata(matdata,i,j)
      if i == j # Diagonal blocks
        _vecdata = select_block_vecdata(vecdata,j)
        _A, _b = A[Block(i,j)], b[Block(i)]
        assemble_matrix_and_vector_add!(_A,_b,_a,(_matvecdata,_matdata,_vecdata))
      else      # Off-diagonal blocks
        __matdata = combine_matdata(_matvecdata,_matdata)
        _A = A[Block(i,j)]
        assemble_matrix_add!(_A,_a,__matdata)
      end
    end
  end
end

function FESpaces.assemble_matrix_and_vector(a::BlockSparseMatrixAssembler,data)
  matvecdata, matdata, vecdata = data
  m = allocate_block_matrix(a)
  v = allocate_block_vector(a)
  block_assemblers = a.block_assemblers
  for i in 1:blocksize(m,1)
    for j in 1:blocksize(m,2)
      _a = block_assemblers[i,j]
      _matvecdata = select_block_matvecdata(matvecdata,i,j)
      _matdata = select_block_matdata(matdata,i,j)
      if i == j # Diagonal blocks
        _vecdata = select_block_vecdata(vecdata,j)
        m[Block(i,j)], v[Block(i)] = assemble_matrix_and_vector(_a,(_matvecdata,_matdata,_vecdata))
      else      # Off-diagonal blocks
        __matdata = combine_matdata(_matvecdata,_matdata) 
        m[Block(i,j)] = assemble_matrix(_a,__matdata)
      end
    end
  end
  return m, v
end
