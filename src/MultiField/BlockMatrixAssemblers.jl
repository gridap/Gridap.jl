
struct BlockMatrixAssembler{A <: FESpaces.Assembler} <: FESpaces.Assembler
  block_assemblers :: AbstractMatrix{A}
end

function FESpaces.get_rows(a::BlockMatrixAssembler)
  row_assemblers = a.block_assemblers[:,1]
  return blockedrange(map(a->length(get_rows(a)),row_assemblers))
end

function FESpaces.get_cols(a::BlockMatrixAssembler)
  col_assemblers = a.block_assemblers[1,:]
  return blockedrange(map(a->length(get_cols(a)),col_assemblers))
end

function FESpaces.get_assembly_strategy(a::BlockMatrixAssembler)
  return get_assembly_strategy(first(a.block_assemblers))
end

function allocate_block_vector(a::BlockMatrixAssembler)
  vector_type = get_vector_type(first(a.block_assemblers))
  vector_lengths = blocklengths(get_rows(a))
  BlockArray(undef_blocks,vector_type,vector_lengths)
end

function allocate_block_matrix(a::BlockMatrixAssembler)
  matrix_type = get_matrix_type(first(a.block_assemblers))
  row_lengths = blocklengths(get_rows(a))
  col_lengths = blocklengths(get_cols(a))
  BlockArray(undef_blocks,matrix_type,row_lengths,col_lengths)
end

# BlockMatrixAssembler for sparse matrices

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

  return BlockMatrixAssembler(block_assemblers)
end

function FESpaces.SparseMatrixAssembler(mat,
                                        vec,
                                        trial::MultiFieldFESpace{<:BlockMultiFieldStyle},
                                        test::MultiFieldFESpace{<:BlockMultiFieldStyle},
                                        strategy::AssemblyStrategy=DefaultAssemblyStrategy())
  return BlockSparseMatrixAssembler(trial,test,SparseMatrixBuilder(mat),ArrayBuilder(vec),strategy)
end

# Block extraction functions

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

function FESpaces.assemble_vector(a::BlockMatrixAssembler,vecdata)
  v = allocate_block_vector(a)
  block_assemblers = a.block_assemblers
  for j in 1:blocksize(v,1)
    _a = block_assemblers[j,1]
    _vecdata = select_block_vecdata(vecdata,j)
    v[Block(j)] = assemble_vector(_a,_vecdata)
  end
  return v
end

function FESpaces.allocate_vector(a::BlockMatrixAssembler,vecdata)
  v = allocate_block_vector(a)
  block_assemblers = a.block_assemblers
  for i in 1:blocksize(v,1)
    _a = block_assemblers[i,1]
    _vecdata = select_block_vecdata(vecdata,i)
    v[Block(i)] = allocate_vector(_a,_vecdata)
  end
  return v
end

function FESpaces.assemble_vector!(b::BlockVector,a::BlockMatrixAssembler,vecdata)
  fill!(b,zero(eltype(b)))
  assemble_vector_add!(b,a,vecdata)
end

function FESpaces.assemble_vector_add!(b::BlockVector,a::BlockMatrixAssembler,vecdata)
  block_assemblers = a.block_assemblers
  for j in 1:blocksize(b,1)
    _a = block_assemblers[j,1]
    _vecdata = select_block_vecdata(vecdata,j)
    bj = b[Block(j)]
    assemble_vector_add!(bj,_a,_vecdata)
  end
end

# Matrix assembly 

function FESpaces.assemble_matrix(a::BlockMatrixAssembler,matdata)
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

function FESpaces.allocate_matrix(a::BlockMatrixAssembler,matdata)
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

function FESpaces.assemble_matrix!(mat::BlockMatrix,a::BlockMatrixAssembler,matdata)
  for (i,j) in blockaxes(mat)
    LinearAlgebra.fillstored!(mat[i,j],zero(eltype(mat)))
  end
  assemble_matrix_add!(mat,a,matdata)
end

function FESpaces.assemble_matrix_add!(mat::BlockMatrix,a::BlockMatrixAssembler,matdata)
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

function FESpaces.allocate_matrix_and_vector(a::BlockMatrixAssembler,data)
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

function FESpaces.assemble_matrix_and_vector!(A::BlockMatrix,b::BlockVector,a::BlockMatrixAssembler,data)
  for (i,j) in blockaxes(A)
    LinearAlgebra.fillstored!(A[i,j],zero(eltype(A)))
  end
  fill!(b,zero(eltype(b)))
  assemble_matrix_and_vector_add!(A,b,a,data)
end

function FESpaces.assemble_matrix_and_vector_add!(A::BlockMatrix,b::BlockVector,a::BlockMatrixAssembler,data)
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

function FESpaces.assemble_matrix_and_vector(a::BlockMatrixAssembler,data)
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
