
function _DiscontinuousFESpace(
  vector_type::Type,
  trian::Triangulation,
  cell_fe::CellFE)

  cell_shapefuns, cell_dof_basis = compute_cell_space(cell_fe,trian)

  cell_dof_ids, nfree = compute_discontinuous_cell_dofs(cell_fe.cell_ctype,cell_fe.ctype_num_dofs)

  ndirichlet = 0
  dirichlet_dof_tag = Int8[]
  dirichlet_cells = Int32[]
  cell_is_dirichlet = Fill(false,num_cells(trian))
  ntags = 0

  metadata = CellConformity(cell_fe)
  UnconstrainedFESpace(
    vector_type,
    nfree,
    ndirichlet,
    cell_dof_ids,
    cell_shapefuns,
    cell_dof_basis,
    cell_is_dirichlet,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags,
    metadata
  )
end

function compute_discontinuous_cell_dofs(cell_to_ctype,ctype_to_nldofs)
  ncells = length(cell_to_ctype)
  ptrs = zeros(Int32,ncells+1)
  for (cell, ctype) in enumerate(cell_to_ctype)
    nldofs = ctype_to_nldofs[ctype]
    ptrs[cell+1] = nldofs
  end
  length_to_ptrs!(ptrs)

  nfree = ptrs[end]-1
  data = collect(Int32,1:nfree)
  
  return Table(data,ptrs), nfree
end

function compute_discontinuous_cell_dofs(
  cell_to_ctype, ctype_to_ldof_to_comp, cell_to_tag, dirichlet_components
)
  ncells = length(cell_to_ctype)

  ptrs = zeros(Int32,ncells+1)
  for (cell, ctype) in enumerate(cell_to_ctype)
    nldofs = length(ctype_to_ldof_to_comp[ctype])
    ptrs[cell+1] = nldofs
  end
  length_to_ptrs!(ptrs)

  ndofs = ptrs[end]-1
  data = zeros(Int32,ndofs)
  dirichlet_dof_tag = zeros(Int8,ndofs)

  nfree = 0
  ndir = 0
  for (cell, ctype) in enumerate(cell_to_ctype)
    ldof_to_comp = ctype_to_ldof_to_comp[ctype]
    tag = cell_to_tag[cell]
    if isequal(tag,UNSET)
      nldofs = length(ldof_to_comp)
      data[ptrs[cell]:ptrs[cell]+nldofs-1] = (nfree+1):(nfree+nldofs)
      ptrs[cell] += nldofs
      nfree += nldofs
    else
      for comp in ldof_to_comp
        if dirichlet_components[comp]
          ndir += 1
          data[ptrs[cell]] = -ndir
          dirichlet_dof_tag[ndir] = tag
        else
          nfree += 1
          data[ptrs[cell]] = nfree
        end
        ptrs[cell] += 1
      end
    end
  end
  rewind_ptrs!(ptrs)

  dirichlet_cells = collect(Int32,findall(!isequal(UNSET),cell_to_tag))
  dirichlet_dof_tag = dirichlet_dof_tag[1:ndir]

  return Table(data,ptrs), nfree, ndir, dirichlet_dof_tag, dirichlet_cells
end
