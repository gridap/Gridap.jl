
function DiscontinuousFESpace(reffes::Vector{<:ReferenceFE}, trian::Triangulation)

  cell_to_ctype = get_cell_type(trian)
  cell_map = get_cell_map(trian)

  cell_dofs, nfree = compute_discontinuous_cell_dofs(reffes,cell_to_ctype)

  ndirichlet = 0
  dirichlet_dof_tag = Int8[]
  dirichlet_cells = Int[]
  ntags = 0

  cell_shapefuns, cell_dof_basis = compute_cell_space(reffes, cell_to_ctype, cell_map)

  UnsconstrainedFESpace(
    nfree,
    ndirichlet,
    cell_dofs,
    cell_shapefuns,
    cell_dof_basis,
    cell_map,
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
  data = collect(Int,1:ndata)

  (Table(data,ptrs), ndata)

end
