
"""
"""
function GradConformingFESpace(reffes,grid,grid_topology,face_labeing,dirichlet_tags)

  T = Float64 # TODO re-think 

  cell_dofs, nfree, ndirichlet, dirichlet_dof_tag = compute_conforming_cell_dofs(
    reffes,grid_topology,face_labeing,dirichlet_tags)

  cell_to_ctype = get_cell_type(grid_topology)
  dof_basis = map(get_dof_basis,reffes)
  cell_dof_basis = CompressedArray(dof_basis,cell_to_ctype)

  shapefuns =  map(get_shapefuns,reffes)
  refshapefuns = CompressedArray(values,cell_to_ctype)
  cell_map = get_cell_map(grid)
  cell_shapefuns = attachmap(refshapefuns,cell_map)

  UnsconstrainedFESpace(
    T,
    nfree,
    ndirichlet,
    cell_dofs,
    cell_shapefuns,
    cell_dof_basis,
    cell_map,
    dirichlet_dof_tag)

end

"""
  compute_conforming_cell_dofs(
    reffes,
    grid_topology,
    face_labeing,
    dirichlet_tags)

  compute_conforming_cell_dofs(
    reffes,
    grid_topology,
    face_labeing,
    dirichlet_tags,
    dirichlet_components)

The result is the tuple

    (cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells)

Assumes that the reffes are aligned with the cell type in the grid_topology
and that it is possible to build a conforming space without imposing constraints

If `dirichlet_components`  is given, then `get_dof_to_comp` has to be defined
for the reference elements in `reffes`.
"""
function compute_conforming_cell_dofs(
  reffes,grid_topology,face_labeing,dirichlet_tags,dirichlet_components=nothing)

  D = num_cell_dims(grid_topology)
  n_faces = num_faces(grid_topology)
  cell_to_ctype = get_cell_type(grid_topology)
  d_to_cell_to_dfaces = [ Table(get_faces(grid_topology,D,d)) for d in 0:D]
  d_to_dface_to_cells = [ Table(get_faces(grid_topology,d,D)) for d in 0:D]
  d_to_offset = get_offsets(grid_topology)
  d_to_ctype_to_ldface_to_own_ldofs = [
    [ get_face_own_dofs(reffe,d) for reffe in reffes] for d in 0:D]

  face_to_own_dofs, ntotal, d_to_dface_to_cell, d_to_dface_to_ldface =  _generate_face_to_own_dofs(
     n_faces,
     cell_to_ctype,
     d_to_cell_to_dfaces,
     d_to_dface_to_cells,
     d_to_offset,
     d_to_ctype_to_ldface_to_own_ldofs)

  d_to_dface_to_tag = [ get_face_tag_index(face_labeing,dirichlet_tags,d)  for d in 0:D]
  cell_to_faces = Table(get_cell_faces(grid_topology))

  nfree, ndiri, diri_dof_tag = _split_face_own_dofs_into_free_and_dirichlet!(
    face_to_own_dofs,
    d_to_offset,
    d_to_dface_to_tag,
    d_to_dface_to_cell,
    d_to_dface_to_ldface,
    cell_to_ctype,
    reffes,
    dirichlet_components)

  cell_to_lface_to_pindex = Table(get_cell_permutations(grid_topology))
  ctype_to_lface_to_own_ldofs = [ get_face_own_dofs(reffe) for reffe in reffes]
  ctype_to_num_dofs = [num_dofs(reffe) for reffe in reffes]
  ctype_to_lface_to_pindex_to_pdofs = [get_face_own_dofs_permutations(reffe) for reffe in reffes]

  cell_dofs = CellDofsNonOriented(
    cell_to_faces,
    cell_to_lface_to_pindex,
    cell_to_ctype,
    ctype_to_lface_to_own_ldofs,
    ctype_to_num_dofs,
    face_to_own_dofs,
    ctype_to_lface_to_pindex_to_pdofs)

  diri_cells = _generate_diri_cells(
    d_to_dface_to_tag,
    d_to_cell_to_dfaces)

  (cell_dofs, nfree, ndiri, diri_dof_tag, diri_cells)
end

function _generate_face_to_own_dofs(
  n_faces,
  cell_to_ctype,
  d_to_cell_to_dfaces::Vector{Table{T,P}},
  d_to_dface_to_cells::Vector{Table{T,P}},
  d_to_offset,
  d_to_ctype_to_ldface_to_own_ldofs) where {T,P}

  face_to_own_dofs_ptrs = zeros(P,n_faces+1)

  D = length(d_to_offset)-1

  icell = 1
  d_to_dface_to_cell = [ get_local_item(d_to_dface_to_cells[d+1],icell)  for d in 0:D ]

  d_to_dface_to_ldface = [
    find_local_index(d_to_dface_to_cell[d+1],d_to_cell_to_dfaces[d+1]) for d in 0:D ]

  for d in 0:D
    cell_to_dfaces = d_to_cell_to_dfaces[d+1]
    dface_to_cells = d_to_dface_to_cells[d+1]
    offset = d_to_offset[d+1]
    ctype_to_ldface_to_own_ldofs = d_to_ctype_to_ldface_to_own_ldofs[d+1]
    ctype_to_ldface_to_num_own_ldofs = map( (x) -> length.(x) ,ctype_to_ldface_to_own_ldofs)
    dface_to_cell_owner = d_to_dface_to_cell[d+1]
    dface_to_ldface = d_to_dface_to_ldface[d+1]

    if any( ctype_to_ldface_to_num_own_ldofs .!= 0)
      _generate_face_to_own_dofs_count_d!(
        face_to_own_dofs_ptrs,
        offset,
        cell_to_ctype,
        dface_to_cell_owner,
        dface_to_ldface,
        ctype_to_ldface_to_num_own_ldofs)
    end
  end

  length_to_ptrs!(face_to_own_dofs_ptrs)

  n_dofs = face_to_own_dofs_ptrs[end]-1
  face_to_own_dofs_data = collect(T(1):T(n_dofs))

  face_to_own_dofs = Table(face_to_own_dofs_data,face_to_own_dofs_ptrs)
  (face_to_own_dofs, n_dofs, d_to_dface_to_cell, d_to_dface_to_ldface)
end

function  _generate_face_to_own_dofs_count_d!(
  face_to_own_dofs_ptrs,
  offset,
  cell_to_ctype,
  dface_to_cell_owner,
  dface_to_ldface,
  ctype_to_ldface_to_num_own_ldofs)

  n_dfaces = length(dface_to_ldface)
  for dface in 1:n_dfaces
    cell = dface_to_cell_owner[dface]
    ldface = dface_to_ldface[dface]
    ctype = cell_to_ctype[cell]
    n_own_ldofs = ctype_to_ldface_to_num_own_ldofs[ctype][ldface]
    face = dface + offset
    face_to_own_dofs_ptrs[face+1] = n_own_ldofs
  end
end

function _split_face_own_dofs_into_free_and_dirichlet!(
  face_to_own_dofs,
  d_to_offset,
  d_to_dface_to_tag,
  d_to_dface_to_cell,
  d_to_dface_to_ldface,
  cell_to_ctype,
  reffes,
  dirichlet_components::Nothing)

  _split_face_own_dofs_into_free_and_dirichlet_generic!(
    face_to_own_dofs,
    d_to_offset,
    d_to_dface_to_tag)

end

function _split_face_own_dofs_into_free_and_dirichlet!(
  face_to_own_dofs,
  d_to_offset,
  d_to_dface_to_tag,
  d_to_dface_to_cell,
  d_to_dface_to_ldface,
  cell_to_ctype,
  reffes,
  dirichlet_components)

  D = length(d_to_dface_to_cell)-1

  d_to_ctype_to_ldface_to_own_ldofs = [
    [ get_face_own_dofs(reffe,d) for reffe in reffes] for d in 0:D]

  ctype_to_ldof_to_comp = [get_dof_to_comp(reffe) for reffe in reffes]

  _split_face_own_dofs_into_free_and_dirichlet_with_components!(
    face_to_own_dofs,
    d_to_offset,
    d_to_dface_to_tag,
    d_to_dface_to_cell,
    d_to_dface_to_ldface,
    cell_to_ctype,
    ctype_to_ldof_to_comp,
    d_to_ctype_to_ldface_to_own_ldofs,
    dirichlet_components)
end

function _split_face_own_dofs_into_free_and_dirichlet_generic!(
  face_to_own_dofs,
  d_to_offset,
  d_to_dface_to_tag)

  D = length(d_to_offset)-1
  nfree = 0
  ndiri = 0
  dirichlet_dof_tag = Int8[]
  for d in 0:D
    dface_to_tag = d_to_dface_to_tag[d+1]
    offset = d_to_offset[d+1]
    for (dface, tag) in enumerate(dface_to_tag)
      face = dface + offset
      pini = face_to_own_dofs.ptrs[face]
      pend = face_to_own_dofs.ptrs[face+1]-1
      if tag == UNSET
        for p in pini:pend
          nfree += 1
          face_to_own_dofs.data[p] = nfree
        end
      else
        for p in pini:pend
          ndiri += 1
          face_to_own_dofs.data[p] = -ndiri
          push!(dirichlet_dof_tag,tag)
        end
      end
    end
  end

  (nfree, ndiri, dirichlet_dof_tag)
end

function _split_face_own_dofs_into_free_and_dirichlet_with_components!(
  face_to_own_dofs,
  d_to_offset,
  d_to_dface_to_tag,
  d_to_dface_to_cell,
  d_to_dface_to_ldface,
  cell_to_ctype,
  ctype_to_ldof_to_comp,
  d_to_ctype_to_ldface_to_own_ldofs,
  dirichlet_components)

  D = length(d_to_offset)-1
  nfree = 0
  ndiri = 0
  dirichlet_dof_tag = Int8[]
  for d in 0:D
    dface_to_tag = d_to_dface_to_tag[d+1]
    offset = d_to_offset[d+1]
    dface_to_cell = d_to_dface_to_cell[d+1]
    dface_to_ldface = d_to_dface_to_ldface[d+1]
    ctype_to_ldface_to_own_ldofs = d_to_ctype_to_ldface_to_own_ldofs[d+1]
    ctype_to_ldface_to_num_own_ldofs = map( (x) -> length.(x) ,ctype_to_ldface_to_own_ldofs)
    if all(ctype_to_ldface_to_num_own_ldofs .== 0)
      continue
    end
    for (dface, diritag) in enumerate(dface_to_tag)
      face = dface + offset
      pini = face_to_own_dofs.ptrs[face]
      pend = face_to_own_dofs.ptrs[face+1]-1
      if diritag == UNSET
        for p in pini:pend
          nfree += 1
          face_to_own_dofs.data[p] = nfree
        end
      else
        cell = dface_to_cell[dface]
        ldface = dface_to_ldface[dface]
        ctype = cell_to_ctype[cell]
        ldof_to_comp = ctype_to_ldof_to_comp[ctype]
        own_ldofs = ctype_to_ldface_to_own_ldofs[ctype][ldface]
        comp_to_isdiri = dirichlet_components[diritag]
        for (fldof, p) in enumerate(pini:pend)
          ldof = own_ldofs[fldof]
          comp = ldof_to_comp[ldof]
          isdiri = comp_to_isdiri[comp]
          if isdiri
            ndiri += 1
            face_to_own_dofs.data[p] = -ndiri
            push!(dirichlet_dof_tag,diritag)
          else
            nfree += 1
            face_to_own_dofs.data[p] = nfree
          end
        end
      end
    end
  end

  (nfree, ndiri, dirichlet_dof_tag)
end

function _generate_diri_cells(
  d_to_dface_to_tag,
  d_to_cell_to_dfaces)

  ncells = length(first(d_to_cell_to_dfaces))
  cell_visited = fill(false,ncells)

  diri_cells = Int[]

  D = length(d_to_dface_to_tag)-1
  for d in 0:D
    dface_to_tag = d_to_dface_to_tag[d+1]
    cell_to_dfaces = Table(d_to_cell_to_dfaces[d+1])
    for cell in 1:length(cell_to_dfaces)
      if cell_visited[cell]
        continue
      end
      pini = cell_to_dfaces.ptrs[cell]
      pend = cell_to_dfaces.ptrs[cell+1]-1
      for p in pini:pend
        dface = cell_to_dfaces.data[p]
        tag = dface_to_tag[dface]
        if tag != UNSET
          push!(diri_cells,cell)
          cell_visited[cell] = true
          break
        end
      end
    end
  end

  diri_cells

end

struct CellDofsNonOriented <:AbstractVector{Vector{Int}}
  cell_to_faces::Table{Int,Int32}
  cell_to_lface_to_pindex::Table{Int8,Int32}
  cell_to_ctype::Vector{Int8}
  ctype_to_lface_to_own_ldofs::Vector{Vector{Vector{Int}}}
  ctype_to_num_dofs::Vector{Int}
  face_to_own_dofs::Table{Int,Int32}
  ctype_to_lface_to_pindex_to_pdofs::Vector{Vector{Vector{Vector{Int}}}}
end

Base.size(a::CellDofsNonOriented) = size(a.cell_to_faces)

Base.IndexStyle(::Type{<:CellDofsNonOriented}) = IndexStyle(Table)

function array_cache(a::CellDofsNonOriented)
  n_dofs = testitem(a.ctype_to_num_dofs)
  T = eltype(eltype(a))
  v = zeros(T,n_dofs)
  CachedArray(v)
end

function getindex!(cache,a::CellDofsNonOriented,cell::Integer)
  ctype = a.cell_to_ctype[cell]
  n_dofs = a.ctype_to_num_dofs[ctype]
  setsize!(cache,(n_dofs,))
  dofs = cache.array
  lface_to_own_ldofs = a.ctype_to_lface_to_own_ldofs[ctype]
  p = a.cell_to_faces.ptrs[cell]-1
  for (lface, own_ldofs) in enumerate(lface_to_own_ldofs)
    face = a.cell_to_faces.data[p+lface]
    pindex = a.cell_to_lface_to_pindex.data[p+lface]
    pdofs = a.ctype_to_lface_to_pindex_to_pdofs[ctype][lface][pindex]

    q = a.face_to_own_dofs.ptrs[face]-1
    for (i,ldof) in enumerate(own_ldofs)
      j = pdofs[i]
      dof = a.face_to_own_dofs.data[q+j]
      dofs[ldof] = dof
    end
  end
  dofs
end

function Base.getindex(a::CellDofsNonOriented,cell::Integer)
  cache = array_cache(a)
  getindex!(cache,a,cell)
end

###struct CellDofsOriented{T,P,V<:AbstractVector} <: AbstractVector{Vector{T}}
###  cell_to_faces::Table{T,P}
###  cell_to_ctype::V
###  ctype_to_lface_to_own_ldofs::Vector{Vector{Vector{Int}}}
###  ctype_to_num_dofs::Vector{Int}
###  face_to_own_dofs::Table{T,P}
###end
###
###Base.size(a::CellDofsOriented) = size(a.cell_to_faces)
###
###Base.IndexStyle(::Type{<:CellDofsOriented}) = IndexStyle(Table)
###
###function array_cache(a::CellDofsOriented)
###  n_dofs = testitem(a.ctype_to_num_dofs)
###  T = eltype(eltype(a))
###  v = zeros(T,n_dofs)
###  CachedArray(v)
###end
###
###function getindex!(cache,a::CellDofsOriented,cell::Integer)
###  ctype = a.cell_to_ctype[cell]
###  n_dofs = a.ctype_to_num_dofs[ctype]
###  setsize!(cache,(n_dofs,))
###  dofs = cache.array
###  lface_to_own_ldofs = a.ctype_to_lface_to_own_ldofs[ctype]
###  p = a.cell_to_faces.ptrs[cell]-1
###  for (lface, own_ldofs) in enumerate(lface_to_own_ldofs)
###    face = a.cell_to_faces.data[p+lface]
###    q = a.face_to_own_dofs.ptrs[face]-1
###    for (i,ldof) in enumerate(own_ldofs)
###      dof = a.face_to_own_dofs.data[q+i]
###      dofs[ldof] = dof
###    end
###  end
###  dofs
###end
###
###function Base.getindex(a::CellDofsOriented,cell::Integer)
###  cache = array_cache(a)
###  getindex!(cache,a,cell)
###end
###

