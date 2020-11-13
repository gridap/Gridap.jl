
function _DiscontinuousFESpace(
  vector_type::Type,
  trian::Triangulation,
  cell_reffe::AbstractArray{<:ReferenceFE},
  cell_shapefuns::CellField,
  cell_dof_basis::CellDof)

  ctype_to_reffe, cell_to_ctype = compress_cell_data(cell_reffe)

  cell_dof_ids, nfree = compute_discontinuous_cell_dofs(ctype_to_reffe,cell_to_ctype)

  ndirichlet = 0
  dirichlet_dof_tag = Int8[]
  dirichlet_cells = Int32[]
  ntags = 0

  UnconstrainedFESpace(
    vector_type,
    nfree,
    ndirichlet,
    cell_dof_ids,
    cell_shapefuns,
    cell_dof_basis,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags)
end

"""
"""
function compute_discontinuous_cell_dofs(reffes,cell_type)

  ctype_to_nldofs = map(num_dofs,reffes)

  _compute_discontinuous_cell_dofs(cell_type,ctype_to_nldofs)

end

function _compute_discontinuous_cell_dofs(cell_to_ctype,ctype_to_nldofs)

  ncells = length(cell_to_ctype)
  ptrs = zeros(Int32,ncells+1)
  for (cell, ctype) in enumerate(cell_to_ctype)
    nldofs = ctype_to_nldofs[ctype]
    ptrs[cell+1] = nldofs
  end

  length_to_ptrs!(ptrs)

  ndata = ptrs[end]-1
  data = collect(Int32,1:ndata)

  (Table(data,ptrs), ndata)

end
