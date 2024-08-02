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
  fill!(b,zero(eltype(b)))
  assemble_vector_add!(b,a,vecdata)
end

function assemble_vector_add!(b,a::SparseMatrixAssembler,vecdata)
  numeric_loop_vector!(b,a,vecdata)
  create_from_nz(b)
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
  LinearAlgebra.fillstored!(mat,zero(eltype(mat)))
  assemble_matrix_add!(mat,a,matdata)
end

function assemble_matrix_add!(mat,a::SparseMatrixAssembler,matdata)
  numeric_loop_matrix!(mat,a,matdata)
  create_from_nz(mat)
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
  m2,v2 = nz_allocation(m1,v1)
  symbolic_loop_matrix_and_vector!(m2,v2,a,data)
  create_from_nz(m2,v2)
end

function assemble_matrix_and_vector!(A,b,a::SparseMatrixAssembler, data)
  LinearAlgebra.fillstored!(A,zero(eltype(A)))
  fill!(b,zero(eltype(b)))
  assemble_matrix_and_vector_add!(A,b,a,data)
end

function assemble_matrix_and_vector_add!(A,b,a::SparseMatrixAssembler,data)
  numeric_loop_matrix_and_vector!(A,b,a,data)
  create_from_nz(A,b)
end

function assemble_matrix_and_vector(a::SparseMatrixAssembler, data)
  m1 = nz_counter(get_matrix_builder(a),(get_rows(a),get_cols(a)))
  v1 = nz_counter(get_vector_builder(a),(get_rows(a),))
  symbolic_loop_matrix_and_vector!(m1,v1,a,data)
  m2,v2 = nz_allocation(m1,v1)
  numeric_loop_matrix_and_vector!(m2,v2,a,data)
  create_from_nz(m2,v2)
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

function symbolic_loop_matrix!(A,a::SparseMatrixAssembler,matdata)
  get_mat(a::Tuple) = a[1]
  get_mat(a) = a
  if LoopStyle(A) == DoNotLoop()
    return A
  end
  strategy = get_assembly_strategy(a)
  for (cellmat,_cellidsrows,_cellidscols) in zip(matdata...)
    cellidsrows = map_cell_rows(strategy,_cellidsrows)
    cellidscols = map_cell_cols(strategy,_cellidscols)
    @assert length(cellidscols) == length(cellidsrows)
    if length(cellidscols) > 0
      rows_cache = array_cache(cellidsrows)
      cols_cache = array_cache(cellidscols)
      mat1 = get_mat(first(cellmat))
      rows1 = getindex!(rows_cache,cellidsrows,1)
      cols1 = getindex!(cols_cache,cellidscols,1)
      touch! = TouchEntriesMap()
      touch_cache = return_cache(touch!,A,mat1,rows1,cols1)
      caches = touch_cache, rows_cache, cols_cache
      _symbolic_loop_matrix!(A,caches,cellidsrows,cellidscols,mat1)
    end
  end
  A
end

@noinline function _symbolic_loop_matrix!(A,caches,cell_rows,cell_cols,mat1)
  touch_cache, rows_cache, cols_cache = caches
  touch! = TouchEntriesMap()
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    evaluate!(touch_cache,touch!,A,mat1,rows,cols)
  end
end

function numeric_loop_matrix!(A,a::SparseMatrixAssembler,matdata)
  strategy = get_assembly_strategy(a)
  for (cellmat,_cellidsrows,_cellidscols) in zip(matdata...)
    cellidsrows = map_cell_rows(strategy,_cellidsrows)
    cellidscols = map_cell_cols(strategy,_cellidscols)
    @assert length(cellidscols) == length(cellidsrows)
    @assert length(cellmat) == length(cellidsrows)
    if length(cellmat) > 0
      rows_cache = array_cache(cellidsrows)
      cols_cache = array_cache(cellidscols)
      vals_cache = array_cache(cellmat)
      mat1 = getindex!(vals_cache,cellmat,1)
      rows1 = getindex!(rows_cache,cellidsrows,1)
      cols1 = getindex!(cols_cache,cellidscols,1)
      add! = AddEntriesMap(+)
      add_cache = return_cache(add!,A,mat1,rows1,cols1)
      caches = add_cache, vals_cache, rows_cache, cols_cache
      _numeric_loop_matrix!(A,caches,cellmat,cellidsrows,cellidscols)
    end
  end
  A
end

@noinline function _numeric_loop_matrix!(mat,caches,cell_vals,cell_rows,cell_cols)
  add_cache, vals_cache, rows_cache, cols_cache = caches
  add! = AddEntriesMap(+)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    evaluate!(add_cache,add!,mat,vals,rows,cols)
  end
end

function symbolic_loop_vector!(b,a::SparseMatrixAssembler,vecdata)
  get_vec(a::Tuple) = a[1]
  get_vec(a) = a
  if LoopStyle(b) == DoNotLoop()
    return b
  end
  strategy = get_assembly_strategy(a)
  for (cellvec,_cellids) in zip(vecdata...)
    cellids = map_cell_rows(strategy,_cellids)
    if length(cellids) > 0
      rows_cache = array_cache(cellids)
      vec1 = get_vec(first(cellvec))
      rows1 = getindex!(rows_cache,cellids,1)
      touch! = TouchEntriesMap()
      touch_cache = return_cache(touch!,b,vec1,rows1)
      caches = touch_cache, rows_cache
      _symbolic_loop_vector!(b,caches,cellids,vec1)
    end
  end
  b
end

@noinline function _symbolic_loop_vector!(A,caches,cellids,vec1)
  touch_cache, rows_cache = caches
  touch! = TouchEntriesMap()
  for cell in 1:length(cellids)
    rows = getindex!(rows_cache,cellids,cell)
    evaluate!(touch_cache,touch!,A,vec1,rows)
  end
end

function numeric_loop_vector!(b,a::SparseMatrixAssembler,vecdata)
  strategy = get_assembly_strategy(a)
  for (cellvec, _cellids) in zip(vecdata...)
    cellids = map_cell_rows(strategy,_cellids)
    if length(cellvec) > 0
      rows_cache = array_cache(cellids)
      vals_cache = array_cache(cellvec)
      vals1 = getindex!(vals_cache,cellvec,1)
      rows1 = getindex!(rows_cache,cellids,1)
      add! = AddEntriesMap(+)
      add_cache = return_cache(add!,b,vals1,rows1)
      caches = add_cache, vals_cache, rows_cache
      _numeric_loop_vector!(b,caches,cellvec,cellids)
    end
  end
  b
end

@noinline function _numeric_loop_vector!(vec,caches,cell_vals,cell_rows)
  add_cache, vals_cache, rows_cache = caches
  @assert length(cell_vals) == length(cell_rows)
  add! = AddEntriesMap(+)
  for cell in 1:length(cell_rows)
    rows = getindex!(rows_cache,cell_rows,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    evaluate!(add_cache,add!,vec,vals,rows)
  end
end

function symbolic_loop_matrix_and_vector!(A,b,a::SparseMatrixAssembler,data)
  if LoopStyle(A) == DoNotLoop()
    return A, b
  end
  matvecdata, matdata, vecdata = data
  strategy = get_assembly_strategy(a)
  for (cellmatvec,_cellidsrows,_cellidscols) in zip(matvecdata...)
    cellidsrows = map_cell_rows(strategy,_cellidsrows)
    cellidscols = map_cell_cols(strategy,_cellidscols)
    @assert length(cellidscols) == length(cellidsrows)
    @assert length(cellmatvec) == length(cellidsrows)
    if length(cellmatvec) > 0
      rows_cache = array_cache(cellidsrows)
      cols_cache = array_cache(cellidscols)
      vals_cache = array_cache(cellmatvec)
      mat1, vec1 = getindex!(vals_cache,cellmatvec,1)
      rows1 = getindex!(rows_cache,cellidsrows,1)
      cols1 = getindex!(cols_cache,cellidscols,1)
      touch! = TouchEntriesMap()
      touch_mat_cache = return_cache(touch!,A,mat1,rows1,cols1)
      touch_vec_cache = return_cache(touch!,b,vec1,rows1)
      caches = touch_mat_cache, touch_vec_cache,rows_cache, cols_cache
      _symbolic_loop_matvec!(A,b,caches,cellidsrows,cellidscols,mat1,vec1)
    end
  end
  symbolic_loop_matrix!(A,a,matdata)
  symbolic_loop_vector!(b,a,vecdata)
  A, b
end 

@noinline function _symbolic_loop_matvec!(A,b,caches,cell_rows,cell_cols,mat1,vec1)
  touch_mat_cache, touch_vec_cache, rows_cache, cols_cache = caches
  touch! = TouchEntriesMap()
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    evaluate!(touch_mat_cache,touch!,A,mat1,rows,cols)
    evaluate!(touch_vec_cache,touch!,b,vec1,rows)
  end
end

function numeric_loop_matrix_and_vector!(A,b,a::SparseMatrixAssembler,data)
  strategy = get_assembly_strategy(a)
  matvecdata, matdata, vecdata = data
  for (cellmatvec,_cellidsrows,_cellidscols) in zip(matvecdata...)
    cellidsrows = map_cell_rows(strategy,_cellidsrows)
    cellidscols = map_cell_cols(strategy,_cellidscols)
    @assert length(cellidscols) == length(cellidsrows)
    @assert length(cellmatvec) == length(cellidsrows)
    if length(cellmatvec) > 0
      rows_cache = array_cache(cellidsrows)
      cols_cache = array_cache(cellidscols)
      vals_cache = array_cache(cellmatvec)
      mat1, vec1 = getindex!(vals_cache,cellmatvec,1)
      rows1 = getindex!(rows_cache,cellidsrows,1)
      cols1 = getindex!(cols_cache,cellidscols,1)
      add! = AddEntriesMap(+)
      add_mat_cache = return_cache(add!,A,mat1,rows1,cols1)
      add_vec_cache = return_cache(add!,b,vec1,rows1)
      caches = add_mat_cache, add_vec_cache, vals_cache, rows_cache, cols_cache
      _numeric_loop_matvec!(A,b,caches,cellmatvec,cellidsrows,cellidscols)
    end
  end
  numeric_loop_matrix!(A,a,matdata)
  numeric_loop_vector!(b,a,vecdata)
  A, b
end

@noinline function _numeric_loop_matvec!(A,b,caches,cell_vals,cell_rows,cell_cols)
  add_mat_cache, add_vec_cache, vals_cache, rows_cache, cols_cache = caches
  add! = AddEntriesMap(+)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    matvals, vecvals = vals
    evaluate!(add_mat_cache,add!,A,matvals,rows,cols)
    evaluate!(add_vec_cache,add!,b,vecvals,rows)
  end
end
