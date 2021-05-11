"""
"""
abstract type SparseMatrixAssembler <: Assembler end

"""
"""
function get_matrix_builder(a::SparseMatrixAssembler)
  @abstractmethod
end

"""
"""
function get_vector_builder(a::SparseMatrixAssembler)
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

get_matrix_type(a::SparseMatrixAssembler) = get_array_type(get_matrix_builder(a))
get_vector_type(a::SparseMatrixAssembler) = get_array_type(get_vector_builder(a))

function allocate_vector(a::SparseMatrixAssembler,vecdata)
  v1 = nz_counter(get_vector_builder(a),(get_rows(a),))
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
  v1 = nz_counter(get_vector_builder(a),(get_rows(a),))
  symbolic_loop_vector!(v1,a,vecdata)
  v2 = nz_allocation(v1)
  numeric_loop_vector!(v2,a,vecdata)
  v3 = create_from_nz(v2)
  v3
end

function allocate_matrix(a::SparseMatrixAssembler,matdata)
  m1 = nz_counter(get_matrix_builder(a),(get_rows(a),get_cols(a)))
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
  m1 = nz_counter(get_matrix_builder(a),(get_rows(a),get_cols(a)))
  symbolic_loop_matrix!(m1,a,matdata)
  m2 = nz_allocation(m1)
  numeric_loop_matrix!(m2,a,matdata)
  m3 = create_from_nz(m2)
  m3
end

function allocate_matrix_and_vector(a::SparseMatrixAssembler,data)
  m1 = nz_counter(get_matrix_builder(a),(get_rows(a),get_cols(a)))
  v1 = nz_counter(get_vector_builder(a),(get_rows(a),))
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
  m1 = nz_counter(get_matrix_builder(a),(get_rows(a),get_cols(a)))
  v1 = nz_counter(get_vector_builder(a),(get_rows(a),))
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
  _ = get_matrix_builder(a)
  _ = get_vector_builder(a)
end

struct GenericSparseMatrixAssembler <: SparseMatrixAssembler
  matrix_builder
  vector_builder
  rows::AbstractUnitRange
  cols::AbstractUnitRange
  strategy::AssemblyStrategy
end

function SparseMatrixAssembler(
  mat,
  vec,
  trial::FESpace,
  test::FESpace,
  strategy::AssemblyStrategy=DefaultAssemblyStrategy())

  rows = get_free_dof_ids(test)
  cols = get_free_dof_ids(trial)
  GenericSparseMatrixAssembler(
    SparseMatrixBuilder(mat),
    ArrayBuilder(vec),
    rows,
    cols,
    strategy)
end

function SparseMatrixAssembler(mat,trial::FESpace,test::FESpace)
  mat_builder = SparseMatrixBuilder(mat)
  T = eltype(get_array_type(mat_builder))
  SparseMatrixAssembler(mat_builder,Vector{T},trial,test)
end

"""
"""
function SparseMatrixAssembler(trial::FESpace,test::FESpace)
  T = get_dof_value_type(trial)
  matrix_type = SparseMatrixCSC{T,Int}
  vector_type = Vector{T}
  SparseMatrixAssembler(matrix_type,vector_type,trial,test)
end

get_rows(a::GenericSparseMatrixAssembler) = a.rows

get_cols(a::GenericSparseMatrixAssembler) = a.cols

get_matrix_builder(a::GenericSparseMatrixAssembler) = a.matrix_builder

get_vector_builder(a::GenericSparseMatrixAssembler) = a.vector_builder

get_assembly_strategy(a::GenericSparseMatrixAssembler) = a.strategy

function symbolic_loop_matrix!(A,a::GenericSparseMatrixAssembler,matdata)
  get_mat(a::Tuple) = a[1]
  get_mat(a) = a
  if LoopStyle(A) == DoNotLoop()
    return A
  end
  for (cellmat,_cellidsrows,_cellidscols) in zip(matdata...)
    cellidsrows = map_cell_rows(a.strategy,_cellidsrows)
    cellidscols = map_cell_cols(a.strategy,_cellidscols)
    rows_cache = array_cache(cellidsrows)
    cols_cache = array_cache(cellidscols)
    @assert length(cellidscols) == length(cellidsrows)
    if length(cellidscols) > 0
      mat1 = get_mat(first(cellmat))
      _symbolic_loop_matrix!(A,rows_cache,cols_cache,cellidsrows,cellidscols,mat1)
    end
  end
  A
end

@noinline function _symbolic_loop_matrix!(A,rows_cache,cols_cache,cell_rows,cell_cols,mat1)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    _symbolic_matrix_entries_at_cell!(A,rows,cols,mat1)
  end
end

@inline function _symbolic_matrix_entries_at_cell!(A,rows,cols,mat1)
  add_entries!(A,nothing,rows,cols)
end

@inline function _symbolic_matrix_entries_at_cell!(
  A,rows::ArrayBlock,cols::ArrayBlock,mat1::ArrayBlock)
  ni,nj = size(mat1.touched)
  for j in 1:nj
    for i in 1:ni
      if mat1.touched[i,j]
        _symbolic_matrix_entries_at_cell!(
          A,rows.array[i],cols.array[j],mat1.array[i,j])
      end
    end
  end
end

function numeric_loop_matrix!(A,a::GenericSparseMatrixAssembler,matdata)
  for (cellmat,_cellidsrows,_cellidscols) in zip(matdata...)
    cellidsrows = map_cell_rows(a.strategy,_cellidsrows)
    cellidscols = map_cell_cols(a.strategy,_cellidscols)
    rows_cache = array_cache(cellidsrows)
    cols_cache = array_cache(cellidscols)
    vals_cache = array_cache(cellmat)
    @assert length(cellidscols) == length(cellidsrows)
    @assert length(cellmat) == length(cellidsrows)
    _numeric_loop_matrix!(
      A,vals_cache,rows_cache,cols_cache,cellmat,cellidsrows,cellidscols)
  end
  A
end

@noinline function _numeric_loop_matrix!(
  mat,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    _numeric_matrix_entries_at_cell!(mat,rows,cols,vals)
  end
end

@inline function _numeric_matrix_entries_at_cell!(mat,rows,cols,vals)
  add_entries!(mat,vals,rows,cols)
end

@inline function _numeric_matrix_entries_at_cell!(
  mat,rows::ArrayBlock,cols::ArrayBlock,vals::ArrayBlock)
  ni, nj = size(vals.touched)
  for j in 1:nj
    for i in 1:ni
      if vals.touched[i,j]
        _numeric_matrix_entries_at_cell!(
          mat,rows.array[i],cols.array[j],vals.array[i,j])
      end
    end
  end
end

function symbolic_loop_vector!(b,a::GenericSparseMatrixAssembler,vecdata)
  @notimplementedif LoopStyle(b) == Loop()
  b
end

function numeric_loop_vector!(b,a::GenericSparseMatrixAssembler,vecdata)
  for (cellvec, _cellids) in zip(vecdata...)
    cellids = map_cell_rows(a.strategy,_cellids)
    rows_cache = array_cache(cellids)
    vals_cache = array_cache(cellvec)
    _numeric_loop_vector!(b,vals_cache,rows_cache,cellvec,cellids)
  end
  b
end

@noinline function _numeric_loop_vector!(vec,vals_cache,rows_cache,cell_vals,cell_rows)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_rows)
    rows = getindex!(rows_cache,cell_rows,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    _numeric_vector_entries_at_cell!(vec,rows,vals)
  end
end

@inline function _numeric_vector_entries_at_cell!(vec,rows,vals)
  add_entries!(vec,vals,rows)
end

@inline function _numeric_vector_entries_at_cell!(
  vec,rows::ArrayBlock,vals::ArrayBlock)
  for i in 1:length(vals.touched)
    if vals.touched[i]
      _numeric_vector_entries_at_cell!(vec,rows.array[i],vals.array[i])
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
  for (cellmatvec,_cellidsrows,_cellidscols) in zip(matvecdata...)
    cellidsrows = map_cell_rows(a.strategy,_cellidsrows)
    cellidscols = map_cell_cols(a.strategy,_cellidscols)
    rows_cache = array_cache(cellidsrows)
    cols_cache = array_cache(cellidscols)
    vals_cache = array_cache(cellmatvec)
    @assert length(cellidscols) == length(cellidsrows)
    @assert length(cellmatvec) == length(cellidsrows)
    _numeric_loop_mavec!(A,b,vals_cache,rows_cache,cols_cache,cellmatvec,cellidsrows,cellidscols)
  end
  numeric_loop_matrix!(A,a,matdata)
  numeric_loop_vector!(b,a,vecdata)
  A, b
end

@noinline function _numeric_loop_mavec!(
  A,b,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    matvals, vecvals = vals
    _numeric_matrix_entries_at_cell!(A,rows,cols,matvals)
    _numeric_vector_entries_at_cell!(b,rows,vecvals)
  end
end

