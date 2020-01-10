
"""
"""
struct ConformingFESpace{T,A,B} <: SingleFieldFESpace
  nfree::Int
  ndirichlet::Int
  cell_dofs::A
  cell_fe::B
  dirichlet_dof_tag::Vector{Int8}

  function ConformingFESpace(
    ::Type{T},
    nfree::Int,
    ndirichlet::Int,
    cell_dofs::AbstractArray,
    cell_fe::SingleFieldCellFE,
    dirichlet_dof_tag::Vector{Int8}) where T

    A = typeof(cell_dofs)
    B = typeof(cell_fe)

    new{T,A,B}(nfree,ndirichlet,cell_dofs,cell_fe,dirichlet_dof_tag)
  end
end


function _compute_conforming_cell_dofs()
end

function num_free_dofs(f::ConformingFESpace)
  f.nfree
end

function num_dirichlet_dofs(f::ConformingFESpace)
  f.ndirichlet
end

function get_cell_dofs(f::ConformingFESpace)
  f.cell_dofs
end

function get_cell_fe(f::ConformingFESpace)
  f.cell_fe
end

function get_dirichlet_dof_tag(f::ConformingFESpace)
  f.dirichlet_dof_tag
end

function zero_free_values(f::ConformingFESpace{T}) where T
  zeros(T,num_free_dofs(f))
end

function zero_dirichlet_values(f::ConformingFESpace{T}) where T
  zeros(T,num_dirichlet_dofs(f))
end

function scatter_free_and_dirichlet_values(f::ConformingFESpace,free_values,dirichlet_values)
  cell_dofs = get_cell_dofs(f)
  LocalToGlobalPosNegArray(cell_dofs,free_vals,dirichlet_values)
end

function gather_free_and_dirichlet_values(f::ConformingFESpace,cell_vals)

  cell_dofs = get_cell_dofs(f)

  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)

  free_vals = zero_free_values(f)
  dirichlet_vals = zero_dirichlet_values(f)

  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs)
end

function  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs)

  for cell in 1:length(cell_vals)
    vals = getindex!(cache_vals,cell_vals)
    dofs = getindex!(cache_dofs,cell_dofs)
    for (i,dof) in enumerate(dofs)
      val = vals[i]
      if dof > 0
        free_vals[dof] = val
      elseif dof < 0
        dirichlet_vals[-dof] = val
      else
        @unreachable "dof ids either positibe or negative, not zero"
      end
    end
  end

end





###function _generate_face_to_own_nodes(model,reffes,conf::RegularConformity)
###
###  D = num_cell_dims(model)
###  cell_to_ctype = get_cell_type(model)
###
###  d_to_ctype_to_ldface_to_own_lnodes = [
###    [ get_face_own_nodes(reffe,d) for reffe in reffes] for d in 0:D]
###
###  d_to_dface_to_cells = [ get_faces(model,d,D) for d in 0:D]
###  d_to_cell_to_dfaces = [ get_faces(model,D,d) for d in 0:D]
###
###  d_to_offset = get_offsets(model)
###
###  n_faces = num_faces(model)
###
###  face_to_own_nodes, n_nodes =  _generate_face_to_own_dofs(
###    n_faces,
###    cell_to_ctype,
###    d_to_cell_to_dfaces,
###    d_to_dface_to_cells,
###    d_to_offset,
###    d_to_ctype_to_ldface_to_own_lnodes)
###
###  (face_to_own_nodes, n_nodes)
###
###end
###
###function _generate_face_to_own_dofs(
###  n_faces,
###  cell_to_ctype,
###  d_to_cell_to_dfaces::Vector{Table{T,P}},
###  d_to_dface_to_cells::Vector{Table{T,P}},
###  d_to_offset,
###  d_to_ctype_to_ldface_to_own_ldofs) where {T,P}
###
###  face_to_own_dofs_ptrs = zeros(P,n_faces+1)
###
###  D = length(d_to_offset)-1
###
###  for d in 0:D
###    cell_to_dfaces = d_to_cell_to_dfaces[d+1]
###    dface_to_cells = d_to_dface_to_cells[d+1]
###    offset = d_to_offset[d+1]
###    ctype_to_ldface_to_own_ldofs = d_to_ctype_to_ldface_to_own_ldofs[d+1]
###    ctype_to_ldface_to_num_own_ldofs = map( (x) -> length.(x) ,ctype_to_ldface_to_own_ldofs)
###    icell = 1
###    dface_to_cell_owner = get_local_item(dface_to_cells,icell)
###    dface_to_ldface = find_local_index(dface_to_cell_owner,cell_to_dfaces)
###
###    if any( ctype_to_ldface_to_num_own_ldofs .!= 0)
###      _generate_face_to_own_dofs_count_d!(
###        face_to_own_dofs_ptrs,
###        offset,
###        cell_to_ctype,
###        dface_to_cell_owner,
###        dface_to_ldface,
###        ctype_to_ldface_to_num_own_ldofs)
###    end
###  end
###
###  length_to_ptrs!(face_to_own_dofs_ptrs)
###  n_dofs = face_to_own_dofs_ptrs[end]-1
###  face_to_own_dofs_data = collect(T(1):T(n_dofs))
###
###  face_to_own_dofs = Table(face_to_own_dofs_data,face_to_own_dofs_ptrs)
###  (face_to_own_dofs, n_dofs)
###end
###
###function  _generate_face_to_own_dofs_count_d!(
###  face_to_own_dofs_ptrs,
###  offset,
###  cell_to_ctype,
###  dface_to_cell_owner,
###  dface_to_ldface,
###  ctype_to_ldface_to_num_own_ldofs)
###
###  n_dfaces = length(dface_to_ldface)
###  for dface in 1:n_dfaces
###    cell = dface_to_cell_owner[dface]
###    ldface = dface_to_ldface[dface]
###    ctype = cell_to_ctype[cell]
###    n_own_ldofs = ctype_to_ldface_to_num_own_ldofs[ctype][ldface]
###    face = dface + offset
###    face_to_own_dofs_ptrs[face+1] = n_own_ldofs
###  end
###end
###
###function _generate_cell_to_nodes(
###  face_to_own_nodes, model, reffes, orientation, conformity)
###  @notimplemented
###end
###
###function _generate_cell_to_nodes(
###  face_to_own_nodes, model, reffes, orientation::Val{true}, ::RegularConformity)
###
###  cell_to_faces = get_cell_faces(model)
###  cell_to_ctype = get_cell_type(model)
###  ctype_lface_to_own_lnodes = [get_face_own_nodes(reffe) for reffe in reffes]
###  ctype_to_num_dofs = map(num_dofs,reffes)
###
###  cell_to_nodes = CellDofsOriented(
###    cell_to_faces,
###    cell_to_ctype,
###    ctype_lface_to_own_lnodes,
###    ctype_to_num_dofs,
###    face_to_own_nodes)
###
###  Table(cell_to_nodes)
###
###end
###
###function _generate_cell_to_nodes(
###  face_to_own_nodes, model, reffes, orientation::Val{false}, ::RegularConformity)
###
###  cell_to_faces = get_cell_faces(model)
###  cell_to_lface_to_pindex = get_cell_perm_indices(model)
###  cell_to_ctype = get_cell_type(model)
###  ctype_to_lface_to_own_lnodes = map( get_face_own_nodes , reffes)
###  ctype_to_num_nodes = map( num_nodes , reffes)
###  facereffes, face_to_ftype = extract_face_reffes(model,reffes)
###  ftype_to_pindex_to_pnodes = map(get_own_dofs_permutations,facereffes)
###
###  cell_to_nodes = CellDofsNonOriented(
###    cell_to_faces,
###    cell_to_lface_to_pindex,
###    cell_to_ctype,
###    ctype_to_lface_to_own_lnodes,
###    ctype_to_num_nodes,
###    face_to_own_nodes,
###    face_to_ftype,
###    ftype_to_pindex_to_pnodes)
###
###  Table(cell_to_nodes)
###
###end
###
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
###struct CellDofsNonOriented{T,P,C<:AbstractVector,F<:AbstractVector} <:AbstractVector{Vector{T}}
###  cell_to_faces::Table{T,P}
###  cell_to_lface_to_pindex::Table{T,P}
###  cell_to_ctype::C
###  ctype_to_lface_to_own_ldofs::Vector{Vector{Vector{Int}}}
###  ctype_to_num_dofs::Vector{Int}
###  face_to_own_dofs::Table{T,P}
###  face_to_ftype::F
###  ftype_to_pindex_to_pdofs::Vector{Vector{Vector{Int}}}
###end
###
###Base.size(a::CellDofsNonOriented) = size(a.cell_to_faces)
###
###Base.IndexStyle(::Type{<:CellDofsNonOriented}) = IndexStyle(Table)
###
###function array_cache(a::CellDofsNonOriented)
###  n_dofs = testitem(a.ctype_to_num_dofs)
###  T = eltype(eltype(a))
###  v = zeros(T,n_dofs)
###  CachedArray(v)
###end
###
###function getindex!(cache,a::CellDofsNonOriented,cell::Integer)
###  ctype = a.cell_to_ctype[cell]
###  n_dofs = a.ctype_to_num_dofs[ctype]
###  setsize!(cache,(n_dofs,))
###  dofs = cache.array
###  lface_to_own_ldofs = a.ctype_to_lface_to_own_ldofs[ctype]
###  p = a.cell_to_faces.ptrs[cell]-1
###  for (lface, own_ldofs) in enumerate(lface_to_own_ldofs)
###    face = a.cell_to_faces.data[p+lface]
###    ftype = a.face_to_ftype[face]
###    pindex = a.cell_to_lface_to_pindex.data[p+lface]
###    pdofs = a.ftype_to_pindex_to_pdofs[ftype][pindex]
###
###    q = a.face_to_own_dofs.ptrs[face]-1
###    for (i,ldof) in enumerate(own_ldofs)
###      j = pdofs[i]
###      dof = a.face_to_own_dofs.data[q+j]
###      dofs[ldof] = dof
###    end
###  end
###  dofs
###end
###
###function Base.getindex(a::CellDofsNonOriented,cell::Integer)
###  cache = array_cache(a)
###  getindex!(cache,a,cell)
###end
