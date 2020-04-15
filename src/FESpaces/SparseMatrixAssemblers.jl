
struct SparseMatrixAssembler{E} <: Assembler
  matrix_type::Type{E}
  test::SingleFieldFESpace
  trial::SingleFieldFESpace
end

"""
"""
function SparseMatrixAssembler(test::SingleFieldFESpace,trial::SingleFieldFESpace)
  SparseMatrixAssembler(SparseMatrixCSC{Float64,Int},test,trial)
end

get_test(a::SparseMatrixAssembler) = a.test

get_trial(a::SparseMatrixAssembler) = a.trial

function allocate_vector(a::SparseMatrixAssembler,term_to_cellidsrows)
  zero_free_values(a.test)
end

function assemble_vector!(b,a::SparseMatrixAssembler,term_to_cellvec,term_to_cellidsrows)
  celldofs = get_cell_dofs(a.test)
  fill!(b,zero(eltype(b)))
  for (cellvec, cellids) in zip(term_to_cellvec,term_to_cellidsrows)
    rows = reindex(celldofs,cellids)
    vals = apply_constraints_vector(a.test,cellvec,cellids)
    rows_cache = array_cache(rows)
    vals_cache = array_cache(vals)
    _assemble_vector!(b,vals_cache,rows_cache,vals,rows)
  end
  b
end

function _assemble_vector!(vec,vals_cache,rows_cache,cell_vals,cell_rows)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_rows)
    rows = getindex!(rows_cache,cell_rows,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    for (i,gid) in enumerate(rows)
      if gid > 0
        vec[gid] += vals[i]
      end
    end
  end
end

function allocate_matrix(a::SparseMatrixAssembler,term_to_cellidsrows, term_to_cellidscols)
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)
  n = 0
  for (cellidsrows,cellidscols) in zip(term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    @assert length(cell_cols) == length(cell_rows)
    n += _count_matrix_entries(a.matrix_type,rows_cache,cols_cache,cell_rows,cell_cols)
  end
  I, J, V = allocate_coo_vectors(a.matrix_type,n)
  nini = 0
  for (cellidsrows,cellidscols) in zip(term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    nini = _allocate_matrix!(a.matrix_type,nini,I,J,rows_cache,cols_cache,cell_rows,cell_cols)
  end
  num_rows = num_free_dofs(a.test)
  num_cols = num_free_dofs(a.trial)
  finalize_coo!(a.matrix_type,I,J,V,num_rows,num_cols)
  sparse_from_coo(a.matrix_type,I,J,V,num_rows,num_cols)
end

@noinline function _count_matrix_entries(::Type{M},rows_cache,cols_cache,cell_rows,cell_cols) where M
  n = 0
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    for gidcol in cols
      if gidcol > 0
        for gidrow in rows
          if gidrow > 0
            if is_entry_stored(M,gidrow,gidcol)
              n += 1
            end
          end
        end
      end
    end
  end
  n
end

@noinline function _allocate_matrix!(a::Type{M},nini,I,J,rows_cache,cols_cache,cell_rows,cell_cols) where M
  n = nini
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    for gidcol in cols
      if gidcol > 0
        for gidrow in rows
          if gidrow > 0
            if is_entry_stored(M,gidrow,gidcol)
              n += 1
              @inbounds I[n] = gidrow
              @inbounds J[n] = gidcol
            end
          end
        end
      end
    end
  end
  n
end

function assemble_matrix!(
  mat,a::SparseMatrixAssembler, term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)
  z = zero(eltype(mat))
  fill_entries!(mat,z)
  assemble_matrix_add!(mat,a,term_to_cellmat,term_to_cellidsrows,term_to_cellidscols)
end

function assemble_matrix_add!(
  mat,a::SparseMatrixAssembler, term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)
  for (cellmat_rc,cellidsrows,cellidscols) in zip(term_to_cellmat,term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    cellmat_r = apply_constraints_matrix_cols(a.trial,cellmat_rc,cellidscols)
    cellmat = apply_constraints_matrix_rows(a.test,cellmat_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cellmat)
    _assemble_matrix!(mat,vals_cache,rows_cache,cols_cache,cellmat,cell_rows,cell_cols)
  end
  mat
end

function _assemble_matrix!(mat,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols)
  @assert length(cell_cols) == length(cell_rows)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    for (j,gidcol) in enumerate(cols)
      if gidcol > 0
        for (i,gidrow) in enumerate(rows)
          if gidrow > 0
            v = vals[i,j]
            add_entry!(mat,v,gidrow,gidcol)
          end
        end
      end
    end
  end
end

function assemble_matrix(
  a::SparseMatrixAssembler, term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)
  n = 0
  for (cellidsrows,cellidscols) in zip(term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    @assert length(cell_cols) == length(cell_rows)
    n += _count_matrix_entries(a.matrix_type,rows_cache,cols_cache,cell_rows,cell_cols)
  end
  I, J, V = allocate_coo_vectors(a.matrix_type,n)
  nini = 0
  for (cellmat_rc,cellidsrows,cellidscols) in zip(term_to_cellmat,term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    cellmat_r = apply_constraints_matrix_cols(a.trial,cellmat_rc,cellidscols)
    cellmat = apply_constraints_matrix_rows(a.test,cellmat_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cellmat)
    @assert length(cell_cols) == length(cell_rows)
    @assert length(cellmat) == length(cell_rows)
    nini = _assemble_matrix_fill!(a.matrix_type,nini,I,J,V,vals_cache,rows_cache,cols_cache,cellmat,cell_rows,cell_cols)
  end
  num_rows = num_free_dofs(a.test)
  num_cols = num_free_dofs(a.trial)
  finalize_coo!(a.matrix_type,I,J,V,num_rows,num_cols)
  sparse_from_coo(a.matrix_type,I,J,V,num_rows,num_cols)
end

@noinline function _assemble_matrix_fill!(::Type{M},nini,I,J,V,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols) where M
  n = nini
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    for (j,gidcol) in enumerate(cols)
      if gidcol > 0
        for (i,gidrow) in enumerate(rows)
          if gidrow > 0
            if is_entry_stored(M,gidrow,gidcol)
              n += 1
              @inbounds v = vals[i,j]
              @inbounds I[n] = gidrow
              @inbounds J[n] = gidcol
              @inbounds V[n] = v
            end
          end
        end
      end
    end
  end
  n
end

function assemble_matrix_and_vector!(A,b,a::SparseMatrixAssembler, matvecdata, matdata, vecdata)
  z = zero(eltype(A))
  fill_entries!(A,z)
  fill!(b,zero(eltype(b)))
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
    _assemble_matrix_and_vector!(A,b,vals_cache,rows_cache,cols_cache,cellmatvec,cell_rows,cell_cols)
  end

  for (cellmat_rc,cellidsrows,cellidscols) in zip(matdata...)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    cellmat_r = apply_constraints_matrix_cols(a.trial,cellmat_rc,cellidscols)
    cellmat = apply_constraints_matrix_rows(a.test,cellmat_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cellmat)
    _assemble_matrix!(A,vals_cache,rows_cache,cols_cache,cellmat,cell_rows,cell_cols)
  end

  for (cellvec, cellids) in zip(vecdata...)
    rows = reindex(celldofs_rows,cellids)
    vals = apply_constraints_vector(a.test,cellvec,cellids)
    rows_cache = array_cache(rows)
    vals_cache = array_cache(vals)
    _assemble_vector!(b,vals_cache,rows_cache,vals,rows)
  end

  A, b
end

function _assemble_matrix_and_vector!(A,b,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols)
  @assert length(cell_cols) == length(cell_rows)
  @assert length(cell_vals) == length(cell_rows)
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    matvals, vecvals = vals
    for (j,gidcol) in enumerate(cols)
      if gidcol > 0
        for (i,gidrow) in enumerate(rows)
          if gidrow > 0
            v = matvals[i,j]
            add_entry!(A,v,gidrow,gidcol)
          end
        end
      end
    end
    for (i,gidrow) in enumerate(rows)
      if gidrow > 0
        bi = vecvals[i]
        b[gidrow] += bi
      end
    end
  end
end

function assemble_matrix_and_vector( a::SparseMatrixAssembler, matvecdata, matdata, vecdata)
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)

  term_to_cellidsrows, term_to_cellidscols,  =  _rearange_cell_ids(matvecdata,matdata,vecdata)

  n = 0
  for (cellidsrows,cellidscols) in zip(term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    @assert length(cell_cols) == length(cell_rows)
    n += _count_matrix_entries(a.matrix_type,rows_cache,cols_cache,cell_rows,cell_cols)
  end

  I, J, V = allocate_coo_vectors(a.matrix_type,n)
  b = zero_free_values(a.test)
  nini = 0

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
    nini = _assemble_matrix_and_vector_fill!(a.matrix_type,nini,I,J,V,b,vals_cache,rows_cache,cols_cache,cellmatvec,cell_rows,cell_cols)
  end

  for (cellmat_rc,cellidsrows,cellidscols) in zip(matdata...)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    cellmat_r = apply_constraints_matrix_cols(a.trial,cellmat_rc,cellidscols)
    cellmat = apply_constraints_matrix_rows(a.test,cellmat_r,cellidsrows)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    vals_cache = array_cache(cellmat)
    @assert length(cell_cols) == length(cell_rows)
    @assert length(cellmat) == length(cell_rows)
    nini = _assemble_matrix_fill!(a.matrix_type,nini,I,J,V,vals_cache,rows_cache,cols_cache,cellmat,cell_rows,cell_cols)
  end

  for (cellvec, cellids) in zip(vecdata...)
    rows = reindex(celldofs_rows,cellids)
    vals = apply_constraints_vector(a.test,cellvec,cellids)
    rows_cache = array_cache(rows)
    vals_cache = array_cache(vals)
    _assemble_vector!(b,vals_cache,rows_cache,vals,rows)
  end

  num_rows = num_free_dofs(a.test)
  num_cols = num_free_dofs(a.trial)
  finalize_coo!(a.matrix_type,I,J,V,num_rows,num_cols)
  A = sparse_from_coo(a.matrix_type,I,J,V,num_rows,num_cols)

  (A, b)
end

@noinline function _assemble_matrix_and_vector_fill!(::Type{M},nini,I,J,V,b,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols) where M
  n = nini
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    vals = getindex!(vals_cache,cell_vals,cell)
    matvals, vecvals = vals
    for (j,gidcol) in enumerate(cols)
      if gidcol > 0
        for (i,gidrow) in enumerate(rows)
          if gidrow > 0
            if is_entry_stored(M,gidrow,gidcol)
              n += 1
              @inbounds v = matvals[i,j]
              @inbounds I[n] = gidrow
              @inbounds J[n] = gidcol
              @inbounds V[n] = v
            end
          end
        end
      end
    end
    for (i,gidrow) in enumerate(rows)
      if gidrow > 0
        bi = vecvals[i]
        b[gidrow] += bi
      end
    end
  end
  n
end
