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

"""
"""
function symbolic_loop_matrix!(A,a::SparseMatrixAssembler,matdata)
  @abstractmethod
end

"""
"""
function symbolic_loop_vector!(b,a::SparseMatrixAssembler,vecdata)
  @abstractmethod
end

"""
"""
function symbolic_loop_matrix_and_vector!(A,b,a::SparseMatrixAssembler,data)
  @abstractmethod
end

"""
"""
function numeric_loop_matrix!(A,a::SparseMatrixAssembler,matdata)
  @abstractmethod
end

"""
"""
function numeric_loop_vector!(b,a::SparseMatrixAssembler,vecdata)
  @abstractmethod
end

"""
"""
function numeric_loop_matrix_and_vector!(A,b,a::SparseMatrixAssembler,data)
  @abstractmethod
end

#"""
#"""
#function count_matrix_nnz_coo(a::SparseMatrixAssembler,matdata)
#  @abstractmethod
#end
#
#"""
#"""
#function count_matrix_and_vector_nnz_coo(a::SparseMatrixAssembler,data)
#  @abstractmethod
#end
#
#"""
#"""
#function fill_matrix_coo_symbolic!(I,J,a::SparseMatrixAssembler,matdata,n=0)
#  @abstractmethod
#end
#
#function fill_matrix_and_vector_coo_symbolic!(I,J,a::SparseMatrixAssembler,data,n=0)
#  @abstractmethod
#end
#
#"""
#"""
#function fill_matrix_coo_numeric!(I,J,V,a::SparseMatrixAssembler,matdata,n=0)
#  @abstractmethod
#end
#
#"""
#"""
#function fill_matrix_and_vector_coo_numeric!(I,J,V,b,a::SparseMatrixAssembler,data,n=0)
#  @abstractmethod
#end

function allocate_vector(a::SparseMatrixAssembler,vecdata)
  v1 = nz_counter(get_vector_type(a),(get_rows(a),))
  symbolic_loop_vector!(v1,a,vecdata)
  v2 = nz_allocation(v1)
  symbolic_loop_vector!(v2,a,vecdata)
  v3 = create_from_nz(v2)
  v3
end

function assemble_vector!(b,a::SparseMatrixAssembler,vecdata)
  fill_entries!(b,zero(eltype(b)))
  assemble_vector_add!(b,a,vecdata)
end

function assemble_vector_add!(b,a::SparseMatrixAssembler,vecdata)
  numeric_loop_vector!(b,a,vecdata)
  b
end

function assemble_vector(a::SparseMatrixAssembler,vecdata)
  v1 = nz_counter(get_vector_type(a),(get_rows(a),))
  symbolic_loop_vector!(v1,a,vecdata)
  v2 = nz_allocation(v1)
  numeric_loop_vector!(v2,a,vecdata)
  v3 = create_from_nz(v2)
  v3
end

function allocate_matrix(a::SparseMatrixAssembler,matdata)
  m1 = nz_counter(get_matrix_type(a),(get_rows(a),get_cols(a)))
  symbolic_loop_matrix!(m1,a,matdata)
  m2 = nz_allocation(m1)
  symbolic_loop_matrix!(m2,a,matdata)
  m3 = create_from_nz(m2)
  m3
end

function assemble_matrix!(mat,a::SparseMatrixAssembler,matdata)
  fill_entries!(mat,zero(eltype(mat)))
  assemble_matrix_add!(mat,a,matdata)
end

function assemble_matrix_add!(mat,a::SparseMatrixAssembler,matdata)
  numeric_loop_matrix!(mat,a,matdata)
  mat
end

function assemble_matrix(a::SparseMatrixAssembler,matdata)
  m1 = nz_counter(get_matrix_type(a),(get_rows(a),get_cols(a)))
  symbolic_loop_matrix!(m1,a,matdata)
  m2 = nz_allocation(m1)
  numeric_loop_matrix!(m2,a,matdata)
  m3 = create_from_nz(m2)
  m3
end

function allocate_matrix_and_vector(a::SparseMatrixAssembler,data)
  m1 = nz_counter(get_matrix_type(a),(get_rows(a),get_cols(a)))
  v1 = nz_counter(get_vector_type(a),(get_rows(a),))
  symbolic_loop_matrix_and_vector!(m1,v1,a,data)
  m2 = nz_allocation(m1)
  v2 = nz_allocation(v1)
  symbolic_loop_matrix_and_vector!(m2,v2,a,data)
  m3 = create_from_nz(m2)
  v3 = create_from_nz(v2)
  m3,v3
end

function assemble_matrix_and_vector!(A,b,a::SparseMatrixAssembler, data)
  fill_entries!(A,zero(eltype(A)))
  fill_entries!(b,zero(eltype(b)))
  assemble_matrix_and_vector_add!(A,b,a,data)
  A, b
end

function assemble_matrix_and_vector_add!(A,b,a::SparseMatrixAssembler,data)
  numeric_loop_matrix_and_vector!(A,b,a,data)
  A, b
end

function assemble_matrix_and_vector(a::SparseMatrixAssembler, data)
  m1 = nz_counter(get_matrix_type(a),(get_rows(a),get_cols(a)))
  v1 = nz_counter(get_vector_type(a),(get_rows(a),))
  symbolic_loop_matrix_and_vector!(m1,v1,a,data)
  m2 = nz_allocation(m1)
  v2 = nz_allocation(v1)
  numeric_loop_matrix_and_vector!(m2,v2,a,data)
  m3 = create_from_nz(m2)
  v3 = create_from_nz(v2)
  m3,v3
end

function test_sparse_matrix_assembler(a::SparseMatrixAssembler,matdata,vecdata,data)
  test_assembler(a,matdata,vecdata,data)
  _ = get_matrix_type(a)
  _ = get_vector_type(a)
end

struct GenericSparseMatrixAssembler{M,V} <: SparseMatrixAssembler
  matrix_type::Type{M}
  vector_type::Type{V}
  rows::AbstractUnitRange
  cols::AbstractUnitRange
  strategy::AssemblyStrategy

  function GenericSparseMatrixAssembler(
    matrix_type::Type{M},
    vector_type::Type{V},
    rows::AbstractUnitRange,
    cols::AbstractUnitRange,
    strategy::AssemblyStrategy) where {M,V}
    new{M,V}(matrix_type,vector_type,rows,cols,strategy)
  end
end

function SparseMatrixAssembler(
  mat::Type,vec::Type,trial::FESpace,test::FESpace,strategy::AssemblyStrategy)
  rows = get_free_dof_ids(test)
  cols = get_free_dof_ids(trial)
  GenericSparseMatrixAssembler(mat,vec,rows,cols,strategy)
end

function SparseMatrixAssembler(mat::Type,vec::Type,trial::FESpace,test::FESpace)
  strategy = DefaultAssemblyStrategy()
  rows = get_free_dof_ids(test)
  cols = get_free_dof_ids(trial)
  GenericSparseMatrixAssembler(mat,vec,rows,cols,strategy)
end

function SparseMatrixAssembler(mat::Type,trial::FESpace,test::FESpace)
  strategy = DefaultAssemblyStrategy()
  rows = get_free_dof_ids(test)
  cols = get_free_dof_ids(trial)
  GenericSparseMatrixAssembler(mat,Vector{eltype(mat)},rows,cols,strategy)
end

"""
"""
function SparseMatrixAssembler(trial::FESpace,test::FESpace)
  T = get_dof_value_type(trial)
  matrix_type = SparseMatrixCSC{T,Int}
  vector_type = Vector{T}
  strategy = DefaultAssemblyStrategy()
  rows = get_free_dof_ids(test)
  cols = get_free_dof_ids(trial)
  GenericSparseMatrixAssembler(matrix_type,vector_type,rows,cols,strategy)
end

get_rows(a::GenericSparseMatrixAssembler) = a.rows

get_cols(a::GenericSparseMatrixAssembler) = a.cols

get_matrix_type(a::GenericSparseMatrixAssembler) = a.matrix_type

get_vector_type(a::GenericSparseMatrixAssembler) = a.vector_type

get_assembly_strategy(a::GenericSparseMatrixAssembler) = a.strategy

function symbolic_loop_matrix!(A,a::GenericSparseMatrixAssembler,matdata)
  if LoopStyle(A) == DoNotLoop()
    return A
  end
  for (cellmat,cellidsrows,cellidscols) in zip(matdata...)
    rows_cache = array_cache(cellidsrows)
    cols_cache = array_cache(cellidscols)
    @assert length(cellidscols) == length(cellidsrows)
    if length(cellidscols) > 0
      mat = first(cellmat)
      Is = _get_block_layout(mat)
      _symbolic_loop_matrix!(A,rows_cache,cols_cache,cellidsrows,cellidscols,a.strategy,Is)
    end
  end
  A
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

@noinline function _symbolic_loop_matrix!(A,rows_cache,cols_cache,cell_rows,cell_cols,strategy,Is)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    _symbolic_matrix_entries_at_cell!(A,rows,cols,strategy,Is)
  end
end

@inline function _symbolic_matrix_entries_at_cell!(A,rows,cols,strategy,Is)
  for gidcol in cols
    if gidcol > 0 && col_mask(strategy,gidcol)
      _gidcol = col_map(strategy,gidcol)
      for gidrow in rows
        if gidrow > 0 && row_mask(strategy,gidrow)
          _gidrow = row_map(strategy,gidrow)
          add_entry!(A,nothing,_gidrow,_gidcol)
        end
      end
    end
  end
end

@inline function _symbolic_matrix_entries_at_cell!(
  A,rows::BlockArrayCoo,cols::BlockArrayCoo,strategy,Is::AbstractArray)
  for (I,Is_next) in Is
    i,j = I.n
    _symbolic_matrix_entries_at_cell!(A,rows[Block(i)],cols[Block(j)],strategy,Is_next)
  end
end

function numeric_loop_matrix!(A,a::GenericSparseMatrixAssembler,matdata)
  for (cellmat,cellidsrows,cellidscols) in zip(matdata...)
    rows_cache = array_cache(cellidsrows)
    cols_cache = array_cache(cellidscols)
    vals_cache = array_cache(cellmat)
    @assert length(cellidscols) == length(cellidsrows)
    @assert length(cellmat) == length(cellidsrows)
    _numeric_loop_matrix!(A,vals_cache,rows_cache,cols_cache,cellmat,cellidsrows,cellidscols,a.strategy)
  end
  A
end

@noinline function _numeric_loop_matrix!(mat,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    _numeric_matrix_entries_at_cell!(mat,rows,cols,vals,strategy)
  end
end

@inline function _numeric_matrix_entries_at_cell!(mat,rows,cols,vals,strategy)
  for (j,gidcol) in enumerate(cols)
    if gidcol > 0 && col_mask(strategy,gidcol)
      _gidcol = col_map(strategy,gidcol)
      for (i,gidrow) in enumerate(rows)
        if gidrow > 0 && row_mask(strategy,gidrow)
          _gidrow = row_map(strategy,gidrow)
          v = vals[i,j]
          add_entry!(mat,v,_gidrow,_gidcol) # TODO reduce the granularity here
        end
      end
    end
  end
end

@inline function _numeric_matrix_entries_at_cell!(
  mat,rows::BlockArrayCoo,cols::BlockArrayCoo,vals::BlockArrayCoo,strategy)
  for I in eachblockid(vals)
    if is_nonzero_block(vals,I)
      i,j = I.n
      _numeric_matrix_entries_at_cell!(mat,rows[Block(i)],cols[Block(j)],vals[I],strategy)
    end
  end
end

function symbolic_loop_vector!(b,a::GenericSparseMatrixAssembler,vecdata)
  @notimplementedif LoopStyle(b) == Loop()
  b
end

function numeric_loop_vector!(b,a::GenericSparseMatrixAssembler,vecdata)
  for (cellvec, cellids) in zip(vecdata...)
    rows_cache = array_cache(cellids)
    vals_cache = array_cache(cellvec)
    _numeric_loop_vector!(b,vals_cache,rows_cache,cellvec,cellids,a.strategy)
  end
  b
end

@noinline function _numeric_loop_vector!(vec,vals_cache,rows_cache,cell_vals,cell_rows,strategy)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_rows)
    rows = getindex!(rows_cache,cell_rows,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    _numeric_vector_entries_at_cell!(vec,rows,vals,strategy)
  end
end

@inline function _numeric_vector_entries_at_cell!(vec,rows,vals,strategy)
  for (i,gid) in enumerate(rows)
    if gid > 0 && row_mask(strategy,gid)
      _gid = row_map(strategy,gid)
      add_entry!(vec,vals[i],_gid)
    end
  end
end

@inline function _numeric_vector_entries_at_cell!(vec,rows::BlockArrayCoo,vals::BlockArrayCoo,strategy)
  for I in eachblockid(vals)
    if is_nonzero_block(vals,I)
      _numeric_vector_entries_at_cell!(vec,rows[I],vals[I],strategy)
    end
  end
end

function symbolic_loop_matrix_and_vector!(A,b,a::GenericSparseMatrixAssembler,data)
  matvecdata, matdata, vecdata = data
  symbolic_loop_matrix!(A,a,matvecdata)
  symbolic_loop_matrix!(A,a,matdata)
  symbolic_loop_vector!(b,a,vecdata)
  A, b
end

function numeric_loop_matrix_and_vector!(A,b,a::GenericSparseMatrixAssembler,data)
  matvecdata, matdata, vecdata = data
  for (cellmatvec,cellidsrows,cellidscols) in zip(matvecdata...)
    rows_cache = array_cache(cellidsrows)
    cols_cache = array_cache(cellidscols)
    vals_cache = array_cache(cellmatvec)
    @assert length(cellidscols) == length(cellidsrows)
    @assert length(cellmatvec) == length(cellidsrows)
    _numeric_loop_mavec!(A,b,vals_cache,rows_cache,cols_cache,cellmatvec,cellidsrows,cellidscols,a.strategy)
  end
  numeric_loop_matrix!(A,a,matdata)
  numeric_loop_vector!(b,a,vecdata)
  A, b
end

@noinline function _numeric_loop_mavec!(A,b,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    matvals, vecvals = vals
    _numeric_matrix_entries_at_cell!(A,rows,cols,matvals,strategy)
    _numeric_vector_entries_at_cell!(b,rows,vecvals,strategy)
  end
end

