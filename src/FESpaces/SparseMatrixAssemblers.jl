
struct SingleFieldSparseMatrixAssembler{M,V} <: SparseMatrixAssembler
  matrix_type::Type{M}
  vector_type::Type{V}
  trial::SingleFieldFESpace
  test::SingleFieldFESpace
  strategy::AssemblyStrategy
end

function SparseMatrixAssembler(
  mat::Type,vec::Type,trial::SingleFieldFESpace,test::SingleFieldFESpace,strategy::AssemblyStrategy)
  SingleFieldSparseMatrixAssembler(mat,vec,trial,test,strategy)
end

function SparseMatrixAssembler(mat::Type,vec::Type,trial::SingleFieldFESpace,test::SingleFieldFESpace)
  strategy = DefaultAssemblyStrategy()
  SingleFieldSparseMatrixAssembler(mat,vec,trial,test,strategy)
end

function SparseMatrixAssembler(mat::Type,trial::SingleFieldFESpace,test::SingleFieldFESpace)
  strategy = DefaultAssemblyStrategy()
  SingleFieldSparseMatrixAssembler(mat,Vector{Float64},trial,test,strategy)
end

"""
"""
function SparseMatrixAssembler(trial::SingleFieldFESpace,test::SingleFieldFESpace)
  matrix_type = SparseMatrixCSC{Float64,Int}
  vector_type = Vector{Float64}
  strategy = DefaultAssemblyStrategy()
  SingleFieldSparseMatrixAssembler(matrix_type,vector_type,trial,test,strategy)
end

get_test(a::SingleFieldSparseMatrixAssembler) = a.test

get_trial(a::SingleFieldSparseMatrixAssembler) = a.trial

get_matrix_type(a::SingleFieldSparseMatrixAssembler) = a.matrix_type

get_vector_type(a::SingleFieldSparseMatrixAssembler) = a.vector_type

get_assembly_strategy(a::SingleFieldSparseMatrixAssembler) = a.strategy

function assemble_vector_add!(b,a::SingleFieldSparseMatrixAssembler,vecdata)
  celldofs = get_cell_dofs(a.test)
  for (cellvec, cellids) in zip(vecdata...)
    rows = reindex(celldofs,cellids)
    vals = apply_constraints_vector(a.test,cellvec,cellids)
    rows_cache = array_cache(rows)
    vals_cache = array_cache(vals)
    _assemble_vector!(b,vals_cache,rows_cache,vals,rows,a.strategy)
  end
  b
end

function _assemble_vector!(vec,vals_cache,rows_cache,cell_vals,cell_rows,strategy)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_rows)
    rows = getindex!(rows_cache,cell_rows,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    for (i,gid) in enumerate(rows)
      if gid > 0 && row_mask(strategy,gid)
        _gid = row_map(strategy,gid)
        add_entry!(vec,vals[i],_gid)
      end
    end
  end
end

function count_matrix_nnz_coo(a::SingleFieldSparseMatrixAssembler,matdata)
  _,term_to_cellidsrows, term_to_cellidscols = matdata
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)
  n = 0
  for (cellidsrows,cellidscols) in zip(term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    @assert length(cell_cols) == length(cell_rows)
    n += _count_matrix_entries(a.matrix_type,rows_cache,cols_cache,cell_rows,cell_cols,a.strategy)
  end

  n
end

@noinline function _count_matrix_entries(::Type{M},rows_cache,cols_cache,cell_rows,cell_cols,strategy) where M
  n = 0
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
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
  end
  n
end

function count_matrix_and_vector_nnz_coo(a::SingleFieldSparseMatrixAssembler,data)
  matvecdata, matdata, vecdata = data
  n = count_matrix_nnz_coo(a,matvecdata)
  n += count_matrix_nnz_coo(a,matdata)
  n
end

function fill_matrix_coo_symbolic!(I,J,a::SingleFieldSparseMatrixAssembler,matdata,n=0)
  _,term_to_cellidsrows, term_to_cellidscols = matdata
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)
  nini = n
  for (cellidsrows,cellidscols) in zip(term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    nini = _allocate_matrix!(a.matrix_type,nini,I,J,rows_cache,cols_cache,cell_rows,cell_cols,a.strategy)
  end
  nini
end

@noinline function _allocate_matrix!(a::Type{M},nini,I,J,rows_cache,cols_cache,cell_rows,cell_cols,strategy) where M
  n = nini
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
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
  end
  n
end

function fill_matrix_and_vector_coo_symbolic!(I,J,a::SingleFieldSparseMatrixAssembler,data,n=0)
  matvecdata, matdata, vecdata = data
  nini = fill_matrix_coo_symbolic!(I,J,a,matvecdata,n)
  nini = fill_matrix_coo_symbolic!(I,J,a,matdata,nini)
  nini
end

function assemble_matrix_add!(mat,a::SingleFieldSparseMatrixAssembler,matdata)

  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)
  for (cellmat_rc,cellidsrows,cellidscols) in zip(matdata...)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    cellmat_r = apply_constraints_matrix_cols(a.trial,cellmat_rc,cellidscols)
    cellmat = apply_constraints_matrix_rows(a.test,cellmat_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cellmat)
    _assemble_matrix!(mat,vals_cache,rows_cache,cols_cache,cellmat,cell_rows,cell_cols,a.strategy)
  end
  mat
end

function _assemble_matrix!(mat,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy)
  @assert length(cell_cols) == length(cell_rows)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
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
end

function fill_matrix_coo_numeric!(I,J,V,a::SingleFieldSparseMatrixAssembler,matdata,n=0)

  nini = n
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)
  for (cellmat_rc,cellidsrows,cellidscols) in zip(matdata...)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    cellmat_r = apply_constraints_matrix_cols(a.trial,cellmat_rc,cellidscols)
    cell_vals = apply_constraints_matrix_rows(a.test,cellmat_r,cellidsrows)
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
  end
  n
end

function assemble_matrix_and_vector_add!(A,b,a::SingleFieldSparseMatrixAssembler, data)

  matvecdata, matdata, vecdata = data
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)

  for (cellmatvec_rc,cellidsrows,cellidscols) in zip(matvecdata...)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    cellmatvec_r = apply_constraints_matrix_and_vector_cols(a.trial,cellmatvec_rc,cellidscols)
    cellmatvec = apply_constraints_matrix_and_vector_rows(a.test,cellmatvec_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cellmatvec)
    _assemble_matrix_and_vector!(A,b,vals_cache,rows_cache,cols_cache,cellmatvec,cell_rows,cell_cols,a.strategy)
  end
  assemble_matrix_add!(A,a,matdata)
  assemble_vector_add!(b,a,vecdata)
  A, b
end

function _assemble_matrix_and_vector!(A,b,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy)
  @assert length(cell_cols) == length(cell_rows)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    matvals, vecvals = vals
    for (j,gidcol) in enumerate(cols)
      if gidcol > 0 && col_mask(strategy,gidcol)
        _gidcol = col_map(strategy,gidcol)
        for (i,gidrow) in enumerate(rows)
          if gidrow > 0 && row_mask(strategy,gidrow)
            _gidrow = row_map(strategy,gidrow)
            v = matvals[i,j]
            add_entry!(A,v,_gidrow,_gidcol)
          end
        end
      end
    end
    for (i,gidrow) in enumerate(rows)
      if gidrow > 0 && row_mask(strategy,gidrow)
        _gidrow = row_map(strategy,gidrow)
        bi = vecvals[i]
        add_entry!(b,bi,_gidrow)
      end
    end
  end
end

function fill_matrix_and_vector_coo_numeric!(I,J,V,b,a::SingleFieldSparseMatrixAssembler,data,n=0)
   
  matvecdata, matdata, vecdata = data
  nini = n

  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)

  for (cellmatvec_rc,cellidsrows,cellidscols) in zip(matvecdata...)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    cellmatvec_r = apply_constraints_matrix_and_vector_cols(a.trial,cellmatvec_rc,cellidscols)
    cellmatvec = apply_constraints_matrix_and_vector_rows(a.test,cellmatvec_r,cellidsrows)
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
    for (j,gidcol) in enumerate(cols)
      if gidcol > 0 && col_mask(strategy,gidcol)
        _gidcol = col_map(strategy,gidcol)
        for (i,gidrow) in enumerate(rows)
          if gidrow > 0 && row_mask(strategy,gidrow)
            _gidrow = row_map(strategy,gidrow)
            if is_entry_stored(M,gidrow,gidcol)
              n += 1
              @inbounds v = matvals[i,j]
              @inbounds I[n] = _gidrow
              @inbounds J[n] = _gidcol
              @inbounds V[n] = v
            end
          end
        end
      end
    end
    for (i,gidrow) in enumerate(rows)
      if gidrow > 0 && row_mask(strategy,gidrow)
        _gidrow = row_map(strategy,gidrow)
        bi = vecvals[i]
        add_entry!(b,bi,_gidrow)
      end
    end
  end
  n
end

