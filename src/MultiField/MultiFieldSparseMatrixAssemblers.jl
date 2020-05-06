
"""
    struct MultiFieldSparseMatrixAssembler{E} <: SparseMatrixAssembler
      # private fields
    end
"""
struct MultiFieldSparseMatrixAssembler{M,V} <: SparseMatrixAssembler
  matrix_type::Type{M}
  vector_type::Type{V}
  test::MultiFieldFESpace
  trial::MultiFieldFESpace
  strategy::AssemblyStrategy
end

function SparseMatrixAssembler(
  matrix_type::Type{<:AbstractSparseMatrix},
  vector_type::Type{<:AbstractVector},
  test::MultiFieldFESpace,
  trial::MultiFieldFESpace)
  strategy = DefaultAssemblyStrategy()
  MultiFieldSparseMatrixAssembler(matrix_type,vector_type,test,trial,strategy)
end

function SparseMatrixAssembler(
  matrix_type::Type{<:AbstractSparseMatrix},
  test::MultiFieldFESpace,
  trial::MultiFieldFESpace)
  SparseMatrixAssembler(matrix_type,Vector{Float64},test,trial)
end

function SparseMatrixAssembler(
  test::MultiFieldFESpace,
  trial::MultiFieldFESpace)
  matrix_type = SparseMatrixCSC{Float64,Int}
  SparseMatrixAssembler(matrix_type,test,trial)
end

function SparseMatrixAssembler(
  matrix_type::Type{<:AbstractSparseMatrix},
  test::Vector{<:SingleFieldFESpace},
  trial::Vector{<:SingleFieldFESpace})
  _test = MultiFieldFESpace(test)
  _trial = MultiFieldFESpace(trial)
  SparseMatrixAssembler(matrix_type,_test,_trial)
end

function SparseMatrixAssembler(
  test::Vector{<:SingleFieldFESpace},
  trial::Vector{<:SingleFieldFESpace})

  matrix_type = SparseMatrixCSC{Float64,Int}
  _test = MultiFieldFESpace(test)
  _trial = MultiFieldFESpace(trial)
  SparseMatrixAssembler(matrix_type,_test,_trial)
end

get_test(a::MultiFieldSparseMatrixAssembler) = a.test

get_trial(a::MultiFieldSparseMatrixAssembler) = a.trial

get_matrix_type(a::MultiFieldSparseMatrixAssembler) = a.matrix_type

get_vector_type(a::MultiFieldSparseMatrixAssembler) = a.vector_type

get_assembly_strategy(a::MultiFieldSparseMatrixAssembler) = a.strategy

function assemble_vector_add!(b,a::MultiFieldSparseMatrixAssembler,term_to_cellvec,term_to_cellidsrows)
  celldofs = get_cell_dofs(a.test)
  for (cellvec, cellids) in zip(term_to_cellvec,term_to_cellidsrows)
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
    _rows = getindex!(rows_cache,cell_rows,cell)
    _vals = getindex!(vals_cache,cell_vals,cell)
    nblocks = length(_vals.blocks)
    for block in 1:nblocks
      field, = _vals.coordinates[block]
      vals = _vals.blocks[block]
      rows = _rows.blocks[field]
      for (i,gid) in enumerate(rows)
        if gid > 0 && row_mask(strategy,gid)
          _gid = row_map(strategy,gid)
          vec[_gid] += vals[i]
        end
      end
    end
  end
end

function count_matrix_nnz_coo(a::MultiFieldSparseMatrixAssembler,term_to_cellidsrows, term_to_cellidscols)
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
    _rows = getindex!(rows_cache,cell_rows,cell)
    _cols = getindex!(cols_cache,cell_cols,cell)
    nfields_rows = length(_rows.blocks)
    nfields_cols = length(_cols.blocks)
    for field_col in 1:nfields_cols
      cols = _cols.blocks[field_col]
      for field_row in 1:nfields_rows
        rows = _rows.blocks[field_row]
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
    end
  end
  n
end

function fill_matrix_coo_symbolic!(I,J,a::MultiFieldSparseMatrixAssembler,term_to_cellidsrows, term_to_cellidscols)
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)
  nini = 0
  for (cellidsrows,cellidscols) in zip(term_to_cellidsrows,term_to_cellidscols)
    cell_rows = reindex(celldofs_rows,cellidsrows)
    cell_cols = reindex(celldofs_cols,cellidscols)
    rows_cache = array_cache(cell_rows)
    cols_cache = array_cache(cell_cols)
    nini = _allocate_matrix!(a.matrix_type,nini,I,J,rows_cache,cols_cache,cell_rows,cell_cols,a.strategy)
  end
end

@noinline function _allocate_matrix!(::Type{M},nini,I,J,rows_cache,cols_cache,cell_rows,cell_cols,strategy) where M
  n = nini
  for cell in 1:length(cell_cols)
    _rows = getindex!(rows_cache,cell_rows,cell)
    _cols = getindex!(cols_cache,cell_cols,cell)
    nfields_rows = length(_rows.blocks)
    nfields_cols = length(_cols.blocks)
    for field_col in 1:nfields_cols
      @inbounds cols = _cols.blocks[field_col]
      for field_row in 1:nfields_rows
        @inbounds rows = _rows.blocks[field_row]
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
    end
  end
  n
end

function assemble_matrix_add!(
  mat,a::MultiFieldSparseMatrixAssembler, term_to_cellmat, term_to_cellidsrows, term_to_cellidscols)

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
    _assemble_matrix!(mat,vals_cache,rows_cache,cols_cache,cellmat,cell_rows,cell_cols,a.strategy)
  end
  mat
end

function _assemble_matrix!(mat,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy)
  for cell in 1:length(cell_cols)
    _rows = getindex!(rows_cache,cell_rows,cell)
    _cols = getindex!(cols_cache,cell_cols,cell)
    _vals = getindex!(vals_cache,cell_vals,cell)
    nblocks = length(_vals.blocks)
    for block in 1:nblocks
      field_row, field_col = _vals.coordinates[block]
      vals = _vals.blocks[block]
      rows = _rows.blocks[field_row]
      cols = _cols.blocks[field_col]
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
end

function fill_matrix_coo_numeric!(
  I,J,V,a::MultiFieldSparseMatrixAssembler,term_to_cellmat,term_to_cellidsrows, term_to_cellidscols,n=0)

  nini = n
  celldofs_rows = get_cell_dofs(a.test)
  celldofs_cols = get_cell_dofs(a.trial)
  for (cellmat_rc,cellidsrows,cellidscols) in zip(term_to_cellmat,term_to_cellidsrows,term_to_cellidscols)
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

@noinline function _fill_matrix!(::Type{M},nini,I,J,V,rows_cache,cols_cache,vals_cache,cell_rows,cell_cols,cell_vals,strategy) where M
  n = nini
  for cell in 1:length(cell_cols)
    _rows = getindex!(rows_cache,cell_rows,cell)
    _cols = getindex!(cols_cache,cell_cols,cell)
    _vals = getindex!(vals_cache,cell_vals,cell)
    nblocks = length(_vals.blocks)
    for block in 1:nblocks
      field_row, field_col = _vals.coordinates[block]
      vals = _vals.blocks[block]
      rows = _rows.blocks[field_row]
      cols = _cols.blocks[field_col]
      for (j,gidcol) in enumerate(cols)
        if gidcol > 0 && col_mask(strategy,gidcol)
          _gidcol = col_map(strategy,gidcol)
          for (i,gidrow) in enumerate(rows)
            if gidrow > 0 && row_mask(strategy,gidrow)
              _gidrow = row_map(strategy,gidrow)
              if is_entry_stored(M,_gidrow,_gidcol)
                n += 1
                @inbounds v = vals[i,j]
                @inbounds I[n] = _gidrow
                @inbounds J[n] = _gidcol
                @inbounds V[n] = v
              end
            end
          end
        end
      end
    end
  end
  n
end

function assemble_matrix_and_vector_add!(
  A,b,a::MultiFieldSparseMatrixAssembler, matvecdata, matdata, vecdata)

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
  assemble_matrix_add!(A,a,matdata...)
  assemble_vector_add!(b,a,vecdata...)
  A, b
end

@noinline function _assemble_matrix_and_vector!(
  mat,vec,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy)

  for cell in 1:length(cell_cols)

    _rows = getindex!(rows_cache,cell_rows,cell)
    _cols = getindex!(cols_cache,cell_cols,cell)
    _vals = getindex!(vals_cache,cell_vals,cell)
    _valsmat, _valsvec = _vals

    nblocks = length(_valsmat.blocks)
    for block in 1:nblocks
      field_row, field_col = _valsmat.coordinates[block]
      valsmat = _valsmat.blocks[block]
      rows = _rows.blocks[field_row]
      cols = _cols.blocks[field_col]
      for (j,gidcol) in enumerate(cols)
        if gidcol > 0 && col_mask(strategy,gidcol)
          _gidcol = col_map(strategy,gidcol)
          for (i,gidrow) in enumerate(rows)
            if gidrow > 0 && row_mask(strategy,gidrow)
              _gidrow = row_map(strategy,gidrow)
              v = valsmat[i,j]
              add_entry!(mat,v,_gidrow,_gidcol)
            end
          end
        end
      end
    end

    nblocks = length(_valsvec.blocks)
    for block in 1:nblocks
      field, = _valsvec.coordinates[block]
      valsvec = _valsvec.blocks[block]
      rows = _rows.blocks[field]
      for (i,gid) in enumerate(rows)
        if gid > 0 && row_mask(strategy,gid)
          _gid = row_map(strategy,gid)
          vec[_gid] += valsvec[i]
        end
      end
    end

  end
end

function fill_matrix_and_vector_coo_numeric!(
  I,J,V,b,a::MultiFieldSparseMatrixAssembler,matvecdata, matdata, vecdata,n=0)

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

  nini = fill_matrix_coo_numeric!(I,J,V,a,matdata...,nini)
  assemble_vector_add!(b,a,vecdata...)

  nini
end


@noinline function _assemble_matrix_and_vector_fill!(
  ::Type{M},nini,I,J,V,b,vals_cache,rows_cache,cols_cache,cell_vals,cell_rows,cell_cols,strategy) where M
  n = nini
  for cell in 1:length(cell_cols)

    _rows = getindex!(rows_cache,cell_rows,cell)
    _cols = getindex!(cols_cache,cell_cols,cell)
    _vals = getindex!(vals_cache,cell_vals,cell)
    _valsmat, _valsvec = _vals

    nblocks = length(_valsmat.blocks)
    for block in 1:nblocks
      field_row, field_col = _valsmat.coordinates[block]
      valsmat = _valsmat.blocks[block]
      rows = _rows.blocks[field_row]
      cols = _cols.blocks[field_col]
      for (j,gidcol) in enumerate(cols)
        if gidcol > 0 && col_mask(strategy,gidcol)
          _gidcol = col_map(strategy,gidcol)
          for (i,gidrow) in enumerate(rows)
            if gidrow > 0 && row_mask(strategy,gidrow)
              _gidrow = row_map(strategy,gidrow)
              if is_entry_stored(M,_gidrow,_gidcol)
                n += 1
                @inbounds v = valsmat[i,j]
                @inbounds I[n] = _gidrow
                @inbounds J[n] = _gidcol
                @inbounds V[n] = v
              end
            end
          end
        end
      end
    end

    nblocks = length(_valsvec.blocks)
    for block in 1:nblocks
      field, = _valsvec.coordinates[block]
      valsvec = _valsvec.blocks[block]
      rows = _rows.blocks[field]
      for (i,gid) in enumerate(rows)
        if gid > 0 && row_mask(strategy,gid)
          _gid = row_map(strategy,gid)
          b[_gid] += valsvec[i]
        end
      end
    end

  end
  n
end

