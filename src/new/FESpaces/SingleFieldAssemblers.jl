
struct SparseMatrixAssembler{E} <: Assembler
  matrix_type::E
  test::SingleFieldFESpace
  trial::SingleFieldFESpace
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
    vals = apply_constraints(a.test,cellvec,cellids)
    rows_cache = array_cache(rows)
    vals_cache = array_cache(vals)
    _assemble_vector!(b,vals_cache,rows_cache,vals,rows)
  end
  b
end

function _assemble_vector!(vec,vals_cache,rows_cache,cell_vals,cell_rows)
  @assert length(cell_vals) == length(cell_rows)
  for i in 1:length(cell_rows)
    rows = getindex!(rows_cache,cell_rows,i)
    vals = getindex!(vals_cache,cell_vals,i)
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
  I, J, V = create_coo_vectors(a.matrix_type)
  for (cellidsrows,cellidscols) in zip(term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    _allocate_matrix!(a.matrix_type,I,J,V,rows_cache,cols_cache,cell_rows,cell_cols)
  end
  num_rows = num_free_dofs(a.test)
  num_cols = num_free_dofs(a.trial)
  finalize_coo!(a.matrix_type,I,J,V,num_rows,num_cols)
  sparse_from_coo(a.matrix_type,I,J,V,num_rows,num_cols)
end

function _allocate_matrix!(M,I,J,V,rows_cache,cols_cache,cell_rows,cell_cols)
  @assert length(cell_cols) == length(cell_rows)
  z = zero(eltype(V))
  for cell in 1:length(cell_cols)
    rows = getindex!(rows_cache,cell_rows,cell)
    cols = getindex!(cols_cache,cell_cols,cell)
    for gidcol in cols
      if gidcol > 0
        for gidrow in rows
          if gidrow > 0
           push_coo!(M, I, J, K, gidrow, gidcol, z)
          end
        end
      end
    end
  end
end

function assemble_matrix!(
  mat,a::SparseMatrixAssembler, term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)
  z = zero(eltype(mat))
  fill_entries!(mat,z)
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
            add_entry!(mat,v,i,j)
          end
        end
      end
    end
  end
end


