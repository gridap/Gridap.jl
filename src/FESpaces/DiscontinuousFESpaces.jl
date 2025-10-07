
"""
    struct DiscontinuousCellConformity <: CellConformity

A CellConformity implementation for discontinuous FE spaces. It assumes that 
all the local DOFs belong to the cell interior, and allows for an arbitrary number of 
DOFs per cell (no compression).
"""
struct DiscontinuousCellConformity{A,B} <: CellConformity
  cell_to_num_dofs::A
  cell_to_ldof_to_comp::B

  function DiscontinuousCellConformity(
    cell_to_num_dofs::AbstractVector{<:Integer},
    cell_to_ldof_to_comp::AbstractVector{<:AbstractVector{<:Integer}}
  )
    @check length(cell_to_num_dofs) == length(cell_to_ldof_to_comp)
    A, B = typeof(cell_to_num_dofs), typeof(cell_to_ldof_to_comp)
    new{A,B}(cell_to_num_dofs,cell_to_ldof_to_comp)
  end
end

function DiscontinuousCellConformity(cell_to_num_dofs::AbstractVector{<:Integer})
  cell_to_ldof_to_comp = lazy_map(n -> fill(Int8(0),n), cell_to_num_dofs)
  DiscontinuousCellConformity(cell_to_num_dofs,cell_to_ldof_to_comp)
end

function DiscontinuousCellConformity(cell_basis::AbstractArray{<:AbstractArray{<:Field}})
  cell_to_num_dofs = collect(Int16,lazy_map(length,cell_basis))
  DiscontinuousCellConformity(cell_to_num_dofs)
end

Geometry.num_cells(c::DiscontinuousCellConformity) = length(c.cell_to_num_dofs)
Geometry.get_cell_type(c::DiscontinuousCellConformity) = Base.OneTo(num_cells(c))

function _DiscontinuousFESpace(
  vector_type::Type, trian::Triangulation, cell_fe::CellFE
)
  cell_conformity = CellConformity(cell_fe)
  cell_shapefuns, cell_dof_basis = compute_cell_space(cell_fe,trian)
  cell_dof_ids, nfree = compute_discontinuous_cell_dofs(cell_conformity)

  ntags = 0
  ndirichlet = 0
  dirichlet_dof_tag = Int8[]
  dirichlet_cells = Int32[]
  cell_is_dirichlet = Fill(false,num_cells(trian))

  metadata = cell_conformity
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

function compute_discontinuous_cell_dofs(cell_conformity::CompressedCellConformity)
  cell_ctype = cell_conformity.cell_ctype
  ctype_num_dofs = cell_conformity.ctype_num_dofs
  return compute_discontinuous_cell_dofs(cell_ctype,ctype_num_dofs)
end

function compute_discontinuous_cell_dofs(cell_conformity::DiscontinuousCellConformity)
  ctype_num_dofs = cell_conformity.cell_to_num_dofs
  cell_ctype = Base.OneTo(length(ctype_num_dofs))
  return compute_discontinuous_cell_dofs(cell_ctype,ctype_num_dofs)
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
  cell_conformity::CompressedCellConformity, cell_to_tag, dirichlet_components
)
  cell_ctype = get_cell_type(cell_conformity)
  ctype_to_nldofs = cell_conformity.ctype_num_dofs
  ctype_ldof_comp = cell_conformity.ctype_ldof_comp
  return compute_discontinuous_cell_dofs(
    cell_ctype, ctype_to_nldofs, ctype_ldof_comp, cell_to_tag, dirichlet_components
  )
end

function compute_discontinuous_cell_dofs(
  cell_conformity::DiscontinuousCellConformity, cell_to_tag, dirichlet_components
)
  cell_ctype = get_cell_type(cell_conformity)
  ctype_to_nldofs = cell_conformity.cell_to_num_dofs
  ctype_ldof_comp = cell_conformity.cell_to_ldof_to_comp
  return compute_discontinuous_cell_dofs(
    cell_ctype, ctype_to_nldofs, ctype_ldof_comp, cell_to_tag, dirichlet_components
  )
end

function compute_discontinuous_cell_dofs(
  cell_to_ctype, ctype_to_nldofs, ctype_to_ldof_to_comp, cell_to_tag, dirichlet_components
)
  ncells = length(cell_to_ctype)

  ptrs = zeros(Int32,ncells+1)
  for (cell, ctype) in enumerate(cell_to_ctype)
    nldofs = ctype_to_nldofs[ctype]
    ptrs[cell+1] = nldofs
  end
  length_to_ptrs!(ptrs)

  ndofs = ptrs[end]-1
  data = zeros(Int32,ndofs)
  dirichlet_dof_tag = zeros(Int8,ndofs)

  nfree = 0
  ndir = 0
  for (cell, ctype) in enumerate(cell_to_ctype)
    nldofs = ctype_to_nldofs[ctype]
    tag = cell_to_tag[cell]
    if isequal(tag,UNSET)
      data[ptrs[cell]:ptrs[cell]+nldofs-1] = (nfree+1):(nfree+nldofs)
      ptrs[cell] += nldofs
      nfree += nldofs
    elseif isnothing(dirichlet_components)
      data[ptrs[cell]:ptrs[cell]+nldofs-1] .= -(ndir+1):-1:-(ndir+nldofs)
      ptrs[cell] += nldofs
      dirichlet_dof_tag[ndir+1:ndir+nldofs] .= tag
      ndir += nldofs
    else
      ldof_to_comp = ctype_to_ldof_to_comp[ctype]
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

function generate_cell_dof_mask(
  scell_conformity::DiscontinuousCellConformity,
  tcell_to_scell::Vector{<:Integer},
  d_to_tcell_to_tdface,
  d_to_tdface_to_mask;
  reverse::Bool=false
)
  D = length(d_to_tcell_to_tdface) - 1
  scell_ndofs = scell_conformity.cell_to_num_dofs
  
  cache = array_cache(d_to_tcell_to_tdface[D+1])
  tcell_dof_mask = Vector{Vector{Bool}}(undef, length(tcell_to_scell))
  for (tcell,scell) in enumerate(tcell_to_scell)
    if scell < 1
      tcell_dof_mask[tcell] = Bool[]
      continue
    end
    tdfaces = getindex!(cache,d_to_tcell_to_tdface[D+1],tcell)
    mask = any(d_to_tdface_to_mask[D+1][tdface] for tdface in tdfaces)
    dof_mask = fill(xor(mask,reverse),scell_ndofs[scell])
    tcell_dof_mask[tcell] = dof_mask
  end

  return tcell_dof_mask
end

function generate_dof_mask(
  scell_conformity::DiscontinuousCellConformity,
  scell_dof_ids::AbstractVector,
  tcell_to_scell::AbstractVector{<:Integer},
  d_to_tcell_to_tdface,
  d_to_tdface_to_mask,
  n_free_dofs;
  reverse::Bool=false
)
  D = length(d_to_tcell_to_tdface) - 1
  dofs_cache = array_cache(scell_dof_ids)

  dof_to_mask = fill(false,n_free_dofs)
  for (tcell,scell) in enumerate(tcell_to_scell)
    (scell < 1) && continue
    if d_to_tdface_to_mask[D+1][tcell]
      dofs = getindex!(dofs_cache,scell_dof_ids,scell)
      for dof in dofs
        if dof > 0 # Avoid dirichlet dofs
          dof_to_mask[dof] = true
        end
      end
    end
  end

  if reverse
    dof_to_mask .= .!dof_to_mask
  end
  return dof_to_mask
end
