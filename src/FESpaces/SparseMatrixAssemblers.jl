# This is an extended interface that only makes sense for assemblers that build (sequential) sparse matrices
# (e.g. not for matrix free assemblers or for distributed assemblers)

"""
"""
abstract type SparseMatrixAssembler <: Assembler end

"""
"""
function get_matrix_type(a::SparseMatrixAssembler)
  @abstractmethod
end

"""
"""
function get_vector_type(a::SparseMatrixAssembler)
  @abstractmethod
end

function allocate_vector(a::SparseMatrixAssembler,vecdata)
  n = num_free_dofs(get_test(a))
  allocate_vector(get_vector_type(a),n)
end

function assemble_vector!(b,a::SparseMatrixAssembler,vecdata)
  fill_entries!(b,zero(eltype(b)))
  assemble_vector_add!(b,a,vecdata)
end

"""
"""
function count_matrix_nnz_coo(a::SparseMatrixAssembler,matdata)
  @abstractmethod
end

"""
"""
function count_matrix_and_vector_nnz_coo(a::SparseMatrixAssembler,data)
  @abstractmethod
end

"""
"""
function fill_matrix_coo_symbolic!(I,J,a::SparseMatrixAssembler,matdata,n=0)
  @abstractmethod
end

function fill_matrix_and_vector_coo_symbolic!(I,J,a::SparseMatrixAssembler,data,n=0)
  @abstractmethod
end

function allocate_matrix(a::SparseMatrixAssembler,matdata)
  n = count_matrix_nnz_coo(a,matdata)
  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)
  fill_matrix_coo_symbolic!(I,J,a,matdata)
  m = num_free_dofs(get_test(a))
  n = num_free_dofs(get_trial(a))
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  sparse_from_coo(get_matrix_type(a),I,J,V,m,n)
end

function assemble_matrix!(mat,a::SparseMatrixAssembler,matdata)
  z = zero(eltype(mat))
  fill_entries!(mat,z)
  assemble_matrix_add!(mat,a,matdata)
end

"""
"""
function fill_matrix_coo_numeric!(I,J,V,a::SparseMatrixAssembler,matdata,n=0)
  @abstractmethod
end

function assemble_matrix(a::SparseMatrixAssembler,matdata)

  n = count_matrix_nnz_coo(a,matdata)
  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)

  fill_matrix_coo_numeric!(I,J,V,a,matdata)

  m = num_free_dofs(get_test(a))
  n = num_free_dofs(get_trial(a))
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  sparse_from_coo(get_matrix_type(a),I,J,V,m,n)
end

function allocate_matrix_and_vector(a::SparseMatrixAssembler,data)

  n = count_matrix_and_vector_nnz_coo(a,data)

  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)
  fill_matrix_and_vector_coo_symbolic!(I,J,a,data)
  m = num_free_dofs(get_test(a))
  n = num_free_dofs(get_trial(a))
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  A = sparse_from_coo(get_matrix_type(a),I,J,V,m,n)

  b = allocate_vector(get_vector_type(a),m)

  A,b
end

function assemble_matrix_and_vector!(A,b,a::SparseMatrixAssembler, data)
  fill_entries!(A,zero(eltype(A)))
  fill_entries!(b,zero(eltype(b)))
  assemble_matrix_and_vector_add!(A,b,a,data)
  A, b
end

"""
"""
function fill_matrix_and_vector_coo_numeric!(I,J,V,b,a::SparseMatrixAssembler,data,n=0)
  @abstractmethod
end

function assemble_matrix_and_vector(a::SparseMatrixAssembler, data)

  n = count_matrix_and_vector_nnz_coo(a,data)
  I,J,V = allocate_coo_vectors(get_matrix_type(a),n)
  n = num_free_dofs(get_test(a))
  b = allocate_vector(get_vector_type(a),n)

  fill_matrix_and_vector_coo_numeric!(I,J,V,b,a,data)

  m = num_free_dofs(get_test(a))
  n = num_free_dofs(get_trial(a))
  finalize_coo!(get_matrix_type(a),I,J,V,m,n)
  A = sparse_from_coo(get_matrix_type(a),I,J,V,m,n)

  A, b
end

function test_sparse_matrix_assembler(a::SparseMatrixAssembler,matdata,vecdata,data)
  test_assembler(a,matdata,vecdata,data)
  _ = get_matrix_type(a)
  _ = get_vector_type(a)
end

struct GenericSparseMatrixAssembler{M,V} <: SparseMatrixAssembler
  matrix_type::Type{M}
  vector_type::Type{V}
  trial::FESpace
  test::FESpace
  strategy::AssemblyStrategy

  function GenericSparseMatrixAssembler(
    matrix_type::Type{M},
    vector_type::Type{V},
    trial::FESpace,
    test::FESpace,
    strategy::AssemblyStrategy) where {M,V}

    @assert ! isa(test,TrialFESpace) """\n
    It is not allowed to build an Assembler with a test space of type TrialFESpace.
    
    Make sure that you are writing first the trial space and then the test space when
    building an Assembler or a FEOperator.
    """
    new{M,V}(matrix_type,vector_type,trial,test,strategy)
  end
end

function SparseMatrixAssembler(
  mat::Type,vec::Type,trial::FESpace,test::FESpace,strategy::AssemblyStrategy)
  GenericSparseMatrixAssembler(mat,vec,trial,test,strategy)
end

function SparseMatrixAssembler(mat::Type,vec::Type,trial::FESpace,test::FESpace)
  strategy = DefaultAssemblyStrategy()
  GenericSparseMatrixAssembler(mat,vec,trial,test,strategy)
end

function SparseMatrixAssembler(mat::Type,trial::FESpace,test::FESpace)
  strategy = DefaultAssemblyStrategy()
  GenericSparseMatrixAssembler(mat,Vector{eltype(mat)},trial,test,strategy)
end

"""
"""
function SparseMatrixAssembler(trial::FESpace,test::FESpace)
  T = get_dof_value_type(trial)
  matrix_type = SparseMatrixCSC{T,Int}
  vector_type = Vector{T}
  strategy = DefaultAssemblyStrategy()
  GenericSparseMatrixAssembler(matrix_type,vector_type,trial,test,strategy)
end

get_test(a::GenericSparseMatrixAssembler) = a.test

get_trial(a::GenericSparseMatrixAssembler) = a.trial

get_matrix_type(a::GenericSparseMatrixAssembler) = a.matrix_type

get_vector_type(a::GenericSparseMatrixAssembler) = a.vector_type

get_assembly_strategy(a::GenericSparseMatrixAssembler) = a.strategy

function assemble_vector_add!(b,a::GenericSparseMatrixAssembler,vecdata)
  for (cellvec, cellids) in zip(vecdata...)
    rows = get_cell_dof_ids(a.test,cellids)
    vals = attach_constraints_rows(a.test,cellvec,cellids)
    rows_cache = array_cache(rows)
    vals_cache = array_cache(vals)
    _assemble_vector!(b,vals_cache,rows_cache,vals,rows,a.strategy)
  end
  b
end

@noinline function _assemble_vector!(vec,vals_cache,rows_cache,cell_vals,cell_rows,strategy)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_rows)
    rows = getindex!(rows_cache,cell_rows,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    _assemble_vector_at_cell!(vec,rows,vals,strategy)
  end
end

@inline function _assemble_vector_at_cell!(vec,rows,vals,strategy)
  for (i,gid) in enumerate(rows)
    if gid > 0 && row_mask(strategy,gid)
      _gid = row_map(strategy,gid)
      add_entry!(vec,vals[i],_gid)
    end
  end
end

@inline function _assemble_vector_at_cell!(vec,rows::BlockArrayCoo,vals::BlockArrayCoo,strategy)
  for I in eachblockid(vals)
    if is_nonzero_block(vals,I)
      _assemble_vector_at_cell!(vec,rows[I],vals[I],strategy)
    end
  end
end

function count_matrix_nnz_coo(a::GenericSparseMatrixAssembler,matdata)
  n = 0
  for (cellmat_rc,cellidsrows,cellidscols) in zip(matdata...)
    cell_rows = get_cell_dof_ids(a.test,cellidsrows)
    cell_cols = get_cell_dof_ids(a.trial,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    cellmat_r = attach_constraints_cols(a.trial,cellmat_rc,cellidscols)
    cellmat = attach_constraints_rows(a.test,cellmat_r,cellidsrows)
    @assert length(cell_cols) == length(cell_rows)
    if length(cell_cols) > 0
      mat = first(cellmat)
      Is = _get_block_layout(mat)
      n += _count_matrix_entries(a.matrix_type,rows_cache,cols_cache,cell_rows,cell_cols,a.strategy,Is)
    end
  end
  n
end

function _get_block_layout(a::Tuple)
  _get_block_layout(a[1])
end

function _get_block_layout(a::AbstractMatrix)
  nothing
end

function _get_block_layout(a::BlockArrayCoo)
  [(I,_get_block_layout(a[I])) for I in eachblockid(a) if is_nonzero_block(a,I) ]
end

@noinline function _count_matrix_entries(::Type{M},rows_cache,cols_cache,cell_rows,cell_cols,strategy,Is) where M
  n = 0
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    n += _count_matrix_entries_at_cell(M,rows,cols,strategy,Is)
  end
  n
end

@inline function _count_matrix_entries_at_cell(::Type{M},rows,cols,strategy,Is) where M
  n = 0
  for gidcol in cols
    if gidcol > 0 && col_mask(strategy,gidcol)
      _gidcol = col_map(strategy,gidcol)
      for gidrow in rows
        if gidrow > 0 && row_mask(strategy,gidrow)
          _gidrow = row_map(strategy,gidrow)
          if is_entry_stored(M,_gidrow,_gidcol)
            n += 1
          end
        end
      end
    end
  end
  n
end

@inline function _count_matrix_entries_at_cell(
  ::Type{M},rows::BlockArrayCoo,cols::BlockArrayCoo,strategy,Is::AbstractArray) where M
  n = 0
  for (I,Is_next) in Is
    i,j = I.n
    n += _count_matrix_entries_at_cell(M,rows[Block(i)],cols[Block(j)],strategy,Is_next)
  end
  n
end

function count_matrix_and_vector_nnz_coo(a::GenericSparseMatrixAssembler,data)
  matvecdata, matdata, vecdata = data
  n = count_matrix_nnz_coo(a,matvecdata)
  n += count_matrix_nnz_coo(a,matdata)
  n
end

function fill_matrix_coo_symbolic!(I,J,a::GenericSparseMatrixAssembler,matdata,n=0)
  term_to_cellmat,term_to_cellidsrows, term_to_cellidscols = matdata
  nini = n
  for (cellmat_rc,cellidsrows,cellidscols) in zip(term_to_cellmat,term_to_cellidsrows,term_to_cellidscols)
    cell_rows = get_cell_dof_ids(a.test,cellidsrows)
    cell_cols = get_cell_dof_ids(a.trial,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    cellmat_r = attach_constraints_cols(a.trial,cellmat_rc,cellidscols)
    cellmat = attach_constraints_rows(a.test,cellmat_r,cellidsrows)
    @assert length(cell_cols) == length(cell_rows)
    if length(cell_cols) > 0
      mat = first(cellmat)
      Is = _get_block_layout(mat)
      nini = _allocate_matrix!(a.matrix_type,nini,I,J,rows_cache,cols_cache,cell_rows,cell_cols,a.strategy,Is)
    end
  end
  nini
end

@noinline function _allocate_matrix!(a::Type{M},nini,I,J,rows_cache,cols_cache,cell_rows,cell_cols,strategy,Is) where M
  n = nini
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    n = _allocate_matrix_at_cell!(M,n,I,J,rows,cols,strategy,Is)
  end
  n
end

@inline function _allocate_matrix_at_cell!(::Type{M},nini,I,J,rows,cols,strategy,Is) where M
  n = nini
  for gidcol in cols
    if gidcol > 0 && col_mask(strategy,gidcol)
      _gidcol = col_map(strategy,gidcol)
      for gidrow in rows
        if gidrow > 0 && row_mask(strategy,gidrow)
          _gidrow = row_map(strategy,gidrow)
          if is_entry_stored(M,_gidrow,_gidcol)
            n += 1
            @inbounds I[n] = _gidrow
            @inbounds J[n] = _gidcol
          end
        end
      end
    end
  end
  n
end

@inline function _allocate_matrix_at_cell!(
  ::Type{M},nini,I,J,rows::BlockArrayCoo,cols::BlockArrayCoo,strategy,Is::AbstractArray) where M
  n = nini
  for (B,Is_next) in Is
    i,j = B.n
    n = _allocate_matrix_at_cell!(M,n,I,J,rows[Block(i)],cols[Block(j)],strategy,Is_next)
  end
  n
end

function fill_matrix_and_vector_coo_symbolic!(I,J,a::GenericSparseMatrixAssembler,data,n=0)
  matvecdata, matdata, vecdata = data
  nini = fill_matrix_coo_symbolic!(I,J,a,matvecdata,n)
  nini = fill_matrix_coo_symbolic!(I,J,a,matdata,nini)
  nini
end

function assemble_matrix_add!(mat,a::GenericSparseMatrixAssembler,matdata)

  for (cellmat_rc,cellidsrows,cellidscols) in zip(matdata...)
    cell_rows = get_cell_dof_ids(a.test,cellidsrows)
    cell_cols = get_cell_dof_ids(a.trial,cellidscols)
    cellmat_r = attach_constraints_cols(a.trial,cellmat_rc,cellidscols)
    cell_vals = attach_constraints_rows(a.test,cellmat_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cell_vals)
    @assert length(cell_cols) == length(cell_rows)
    @assert length(cell_vals) == length(cell_rows)
    _assemble_matrix!(mat,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,a.strategy)
  end
  mat
end

@noinline function _assemble_matrix!(mat,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    _assemble_matrix_at_cell!(mat,rows,cols,vals,strategy)
  end
end

@inline function _assemble_matrix_at_cell!(mat,rows,cols,vals,strategy)
  for (j,gidcol) in enumerate(cols)
    if gidcol > 0 && col_mask(strategy,gidcol)
      _gidcol = col_map(strategy,gidcol)
      for (i,gidrow) in enumerate(rows)
        if gidrow > 0 && row_mask(strategy,gidrow)
          _gidrow = row_map(strategy,gidrow)
          v = vals[i,j]
          add_entry!(mat,v,_gidrow,_gidcol)
        end
      end
    end
  end
end

@inline function _assemble_matrix_at_cell!(mat,rows::BlockArrayCoo,cols::BlockArrayCoo,vals::BlockArrayCoo,strategy)
  for I in eachblockid(vals)
    if is_nonzero_block(vals,I)
      i,j = I.n
      _assemble_matrix_at_cell!(mat,rows[Block(i)],cols[Block(j)],vals[I],strategy)
    end
  end
end

function fill_matrix_coo_numeric!(I,J,V,a::GenericSparseMatrixAssembler,matdata,n=0)
  nini = n
  for (cellmat_rc,cellidsrows,cellidscols) in zip(matdata...)
    cell_rows = get_cell_dof_ids(a.test,cellidsrows)
    cell_cols = get_cell_dof_ids(a.trial,cellidscols)
    cellmat_r = attach_constraints_cols(a.trial,cellmat_rc,cellidscols)
    cell_vals = attach_constraints_rows(a.test,cellmat_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cell_vals)
    nini = _fill_matrix!(
      a.matrix_type,nini,I,J,V,rows_cache,cols_cache,vals_cache,cell_rows,cell_cols,cell_vals,a.strategy)
  end

  nini
end

@noinline function _fill_matrix!(
  a::Type{M},nini,I,J,V,rows_cache,cols_cache,vals_cache,cell_rows,cell_cols,cell_vals,strategy) where M

  n = nini
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    n = _fill_matrix_at_cell!(M,n,I,J,V,rows,cols,vals,strategy)
  end
  n
end

@inline function _fill_matrix_at_cell!(::Type{M},nini,I,J,V,rows,cols,vals,strategy) where M
  n = nini
  for (j,gidcol) in enumerate(cols)
    if gidcol > 0 && col_mask(strategy,gidcol)
      _gidcol = col_map(strategy,gidcol)
      for (i,gidrow) in enumerate(rows)
        if gidrow > 0 && row_mask(strategy,gidrow)
          _gidrow = row_map(strategy,gidrow)
          if is_entry_stored(M,_gidrow,_gidcol)
            n += 1
            @inbounds I[n] = _gidrow
            @inbounds J[n] = _gidcol
            @inbounds V[n] = vals[i,j]
          end
        end
      end
    end
  end
  n
end

@inline function _fill_matrix_at_cell!(
  ::Type{M},nini,I,J,V,rows::BlockArrayCoo,cols::BlockArrayCoo,vals::BlockArrayCoo,strategy) where M
  n = nini
  for B in eachblockid(vals)
    if is_nonzero_block(vals,B)
      i,j = B.n
      n = _fill_matrix_at_cell!(M,n,I,J,V,rows[Block(i)],cols[Block(j)],vals[B],strategy)
    end
  end
  n
end

function assemble_matrix_and_vector_add!(A,b,a::GenericSparseMatrixAssembler, data)

  matvecdata, matdata, vecdata = data

  for (cellmatvec_rc,cellidsrows,cellidscols) in zip(matvecdata...)
    cell_rows = get_cell_dof_ids(a.test,cellidsrows)
    cell_cols = get_cell_dof_ids(a.trial,cellidscols)
    cellmatvec_r = attach_constraints_cols(a.trial,cellmatvec_rc,cellidscols)
    cellmatvec = attach_constraints_rows(a.test,cellmatvec_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cellmatvec)
    _assemble_matrix_and_vector!(A,b,vals_cache,rows_cache,cols_cache,cellmatvec,cell_rows,cell_cols,a.strategy)
  end
  assemble_matrix_add!(A,a,matdata)
  assemble_vector_add!(b,a,vecdata)
  A, b
end

@noinline function _assemble_matrix_and_vector!(A,b,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy)
  @assert length(cell_cols) == length(cell_rows)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    matvals, vecvals = vals
    _assemble_matrix_at_cell!(A,rows,cols,matvals,strategy)
    _assemble_vector_at_cell!(b,rows,vecvals,strategy)
  end
end

function fill_matrix_and_vector_coo_numeric!(I,J,V,b,a::GenericSparseMatrixAssembler,data,n=0)
   
  matvecdata, matdata, vecdata = data
  nini = n

  for (cellmatvec_rc,cellidsrows,cellidscols) in zip(matvecdata...)
    cell_rows = get_cell_dof_ids(a.test,cellidsrows)
    cell_cols = get_cell_dof_ids(a.trial,cellidscols)
    cellmatvec_r = attach_constraints_cols(a.trial,cellmatvec_rc,cellidscols)
    cellmatvec = attach_constraints_rows(a.test,cellmatvec_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cellmatvec)
    @assert length(cell_cols) == length(cell_rows)
    @assert length(cellmatvec) == length(cell_rows)
    nini = _assemble_matrix_and_vector_fill!(
      a.matrix_type,nini,I,J,V,b,vals_cache,rows_cache,cols_cache,cellmatvec,cell_rows,cell_cols,a.strategy)
  end

  nini = fill_matrix_coo_numeric!(I,J,V,a,matdata,nini)
  assemble_vector_add!(b,a,vecdata)

  nini
end

@noinline function _assemble_matrix_and_vector_fill!(
  ::Type{M},nini,I,J,V,b,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy) where M
  n = nini
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    matvals, vecvals = vals
    n = _fill_matrix_at_cell!(M,n,I,J,V,rows,cols,matvals,strategy)
    _assemble_vector_at_cell!(b,rows,vecvals,strategy)
  end
  n
end

