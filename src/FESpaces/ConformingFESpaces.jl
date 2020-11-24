
"""
Minimum data required to describe dof ownership.
At this moment, the cell-wise ownership is compressed on cell types.
This can be relaxed in the future, to have an arbitrary cell-wise dof ownership.
"""
struct CellConformity{T} <: GridapType
  cell_ctype::T
  ctype_lface_own_ldofs::Vector{Vector{Vector{Int}}}
  ctype_lface_pindex_pdofs::Vector{Vector{Vector{Vector{Int}}}}
  d_ctype_num_dfaces::Vector{Vector{Int}}
end

function Base.getproperty(a::CellConformity, sym::Symbol)
  if sym == :d_ctype_offset
    _d_ctype_offset(a)
  elseif sym == :d_ctype_ldface_own_ldofs
    _d_ctype_ldface_own_ldofs(a)
  else
    getfield(a, sym)
  end
end

function Base.propertynames(x::CellConformity, private=false)
  (fieldnames(typeof(x))...,:d_ctype_offset,:d_ctype_ldface_own_ldofs)
end

function _d_ctype_offset(a::CellConformity)
  num_ctypes = length(a.ctype_lface_own_ldofs)
  num_ds = length(a.d_ctype_num_dfaces)
  d_ctype_offset = [ zeros(Int,num_ctypes) for d in 1:num_ds]
  for d in 2:num_ds
    for ctype in 1:num_ctypes
      d_ctype_offset[d][ctype] = d_ctype_offset[d-1][ctype] + a.d_ctype_num_dfaces[d-1][ctype]
    end
  end
  d_ctype_offset
end

function _d_ctype_ldface_own_ldofs(a::CellConformity)
  num_ctypes = length(a.ctype_lface_own_ldofs)
  num_ds = length(a.d_ctype_num_dfaces)
  [[[ a.ctype_lface_own_ldofs[ctype][ldface+a.d_ctype_offset[d][ctype]]
    for ldface in 1:a.d_ctype_num_dfaces[d][ctype] ]
    for ctype in 1:num_ctypes ]
    for d in 1:num_ds ]
end

"""
Minimum data required to build a conforming FE space.
At this moment, the some cell-wise info is compressed on cell types.
This can be relaxed in the future, and have an arbitrary cell-wise data.
"""
struct CellFE{T} <: GridapType
  cell_ctype::T
  ctype_num_dofs::Vector{Int}
  ctype_ldof_comp::Vector{Vector{Int}}
  cell_conformity::CellConformity{T}
  cell_shapefuns::AbstractArray{<:AbstractVector{<:Field}}
  cell_dof_basis::AbstractArray{<:AbstractVector{<:Dof}}
  cell_shapefuns_domain::DomainStyle
  cell_dof_basis_domain::DomainStyle
  max_order::Int 
end
# If the shapefuns are not polynomials, max_order has to be understood as the order of a
# reasonable quadrature rule to integrate the shape functions. Only used by FESpace 
# constructors that need to integrate the shape functions (e.g., ZeroMeanFESpace).

Geometry.num_cells(cell_fe::CellFE) = length(cell_fe.cell_ctype)

function CellConformity(cell_fe::CellFE)
  cell_fe.cell_conformity
end

function CellConformity(cell_fe::CellFE,cell_conf::Nothing)
  cell_fe.cell_conformity
end

function CellConformity(cell_fe::CellFE,cell_conf::CellConformity)
  @assert length(cell_fe.cell_ctype) == length(cell_fe.cell_ctype)
  cell_conf
end

function CellConformity(cell_fe::CellFE,conformity::Union{Symbol,Conformity})
  if conformity in (L2Conformity(),:L2)
    @notimplemented
  elseif conformity == :default
    cell_fe.cell_conformity
  else
    @unreachable """\n
    Argument conformity should be either :default, :L2, or L2Conformity()
    """
  end
end

"""
Generate A CellConformity from a vector of reference fes
"""
function CellConformity(cell_reffe::AbstractArray{<:ReferenceFE},conformity=nothing)
  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
  conf = Conformity(first(ctype_reffe),conformity)
  ctype_lface_own_ldofs = map(reffe->get_face_own_dofs(reffe,conf),ctype_reffe)
  ctype_lface_pindex_pdofs = map(reffe->get_face_own_dofs_permutations(reffe,conf),ctype_reffe)
  D = num_dims(first(ctype_reffe))
  d_ctype_num_dfaces = [ map(reffe->num_faces(get_polytope(reffe),d),ctype_reffe) for d in 0:D]
  CellConformity(
    cell_ctype,
    ctype_lface_own_ldofs,
    ctype_lface_pindex_pdofs,
    d_ctype_num_dfaces)
end

"""
Generate a CellFE from a vector of reference fes
"""
function CellFE(
  cell_map::AbstractArray{<:Field},
  cell_reffe::AbstractArray{<:ReferenceFE})

  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
  ctype_num_dofs = map(num_dofs,ctype_reffe)
  ctype_ldof_comp = map(reffe->get_dof_to_comp(reffe),ctype_reffe)
  cell_conformity = CellConformity(cell_reffe)
  cell_shapefuns = lazy_map(get_shapefuns,cell_reffe,cell_map)
  cell_dof_basis = lazy_map(get_dof_basis,cell_reffe,cell_map)
  cell_shapefuns_domain = ReferenceDomain()
  cell_dof_basis_domain = cell_shapefuns_domain
  max_order = maximum(map(get_order,ctype_reffe))

  CellFE(
    cell_ctype,
    ctype_num_dofs,
    ctype_ldof_comp,
    cell_conformity,
    cell_shapefuns,
    cell_dof_basis,
    cell_shapefuns_domain,
    cell_dof_basis_domain,
    max_order)
end

# Low level conforming FE Space constructor
# The user is not expected to call this function. Use the factory function instead
function _ConformingFESpace(
  vector_type::Type,
  model::DiscreteModel,
  face_labeling::FaceLabeling,
  cell_fe::CellFE,
  cell_conformity::CellConformity,
  dirichlet_tags,
  dirichlet_components)

  grid_topology = get_grid_topology(model)
  ntags = length(dirichlet_tags)

  cell_dofs_ids, nfree, ndirichlet, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
    cell_fe,cell_conformity,grid_topology,face_labeling,dirichlet_tags,dirichlet_components)

  trian = Triangulation(model)
  cell_shapefuns, cell_dof_basis = compute_cell_space(cell_fe,trian)

  UnconstrainedFESpace(
    vector_type,
    nfree,
    ndirichlet,
    cell_dofs_ids,
    cell_shapefuns,
    cell_dof_basis,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags)
end

function compute_cell_space(cell_fe,trian::Triangulation)
  cell_shapefuns, cell_dof_basis, d1, d2 = _compute_cell_space(cell_fe,trian)
  FEBasis(cell_shapefuns,trian,TestBasis(),d1), CellDof(cell_dof_basis,trian,d2)
end

function _compute_cell_space(cell_fe::CellFE,trian::Triangulation)
  ( cell_fe.cell_shapefuns,
    cell_fe.cell_dof_basis,
    cell_fe.cell_shapefuns_domain,
    cell_fe.cell_dof_basis_domain)
end

function _compute_cell_space(
  cell_reffe::AbstractArray{<:ReferenceFE}, trian::Triangulation)

  cell_map = get_cell_map(trian)
  cell_shapefuns = lazy_map(get_shapefuns,cell_reffe,cell_map)
  cell_dof_basis = lazy_map(get_dof_basis,cell_reffe,cell_map)
  cell_shapefuns, cell_dof_basis, ReferenceDomain(), ReferenceDomain()
end

"""
The result is the tuple

    (cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells)

If `dirichlet_components`  is given, then `get_dof_to_comp` has to be defined
for the reference elements in `reffes`.
"""
function compute_conforming_cell_dofs(
  cell_fe,cell_conformity,grid_topology,face_labeling,dirichlet_tags,dirichlet_components=nothing)

  cell_to_ctype = cell_fe.cell_ctype
  ctype_to_ldof_to_comp = cell_fe.ctype_ldof_comp
  ctype_to_num_dofs = cell_fe.ctype_num_dofs

  d_to_ctype_to_ldface_to_own_ldofs = cell_conformity.d_ctype_ldface_own_ldofs
  ctype_to_lface_to_own_ldofs = cell_conformity.ctype_lface_own_ldofs
  ctype_to_lface_to_pindex_to_pdofs = cell_conformity.ctype_lface_pindex_pdofs

  D = num_cell_dims(grid_topology)
  n_faces = num_faces(grid_topology)
  d_to_cell_to_dfaces = [ Table(get_faces(grid_topology,D,d)) for d in 0:D]
  d_to_dface_to_cells = [ Table(get_faces(grid_topology,d,D)) for d in 0:D]
  d_to_offset = get_offsets(grid_topology)

  face_to_own_dofs, ntotal, d_to_dface_to_cell, d_to_dface_to_ldface =  _generate_face_to_own_dofs(
     n_faces,
     cell_to_ctype,
     d_to_cell_to_dfaces,
     d_to_dface_to_cells,
     d_to_offset,
     d_to_ctype_to_ldface_to_own_ldofs)

  d_to_dface_to_tag = [ get_face_tag_index(face_labeling,dirichlet_tags,d)  for d in 0:D]
  cell_to_faces = Table(get_cell_faces(grid_topology))

  _dirichlet_components = _convert_dirichlet_components(dirichlet_tags,dirichlet_components)

  nfree, ndiri, diri_dof_tag = _split_face_own_dofs_into_free_and_dirichlet!(
    face_to_own_dofs,
    d_to_offset,
    d_to_dface_to_tag,
    d_to_dface_to_cell,
    d_to_dface_to_ldface,
    cell_to_ctype,
    d_to_ctype_to_ldface_to_own_ldofs,
    ctype_to_ldof_to_comp,
    _dirichlet_components)

  cell_to_lface_to_pindex = Table(get_cell_permutations(grid_topology))

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

function _convert_dirichlet_components(dirichlet_tags::AbstractArray,dirichlet_components)
  dirichlet_components
end

function _convert_dirichlet_components(dirichlet_tags::Integer,dirichlet_components)
  dirichlet_components==nothing ? nothing : [dirichlet_components,]
end

function _convert_dirichlet_components(dirichlet_tags::AbstractString,dirichlet_components)
  dirichlet_components==nothing ? nothing : [dirichlet_components,]
end

function _generate_face_to_own_dofs(
  n_faces,
  cell_to_ctype,
  d_to_cell_to_dfaces::Vector{Table{T,Vd,Vp}},
  d_to_dface_to_cells::Vector{Table{T,Vd,Vp}},
  d_to_offset,
  d_to_ctype_to_ldface_to_own_ldofs) where {T,Vd,Vp}

  P=eltype(Vp)
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
  d_to_ctype_to_ldface_to_own_ldofs,
  ctype_to_ldof_to_comp,
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
  d_to_ctype_to_ldface_to_own_ldofs,
  ctype_to_ldof_to_comp,
  dirichlet_components)

  D = length(d_to_dface_to_cell)-1

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

  diri_cells = Int32[]

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

struct CellDofsNonOriented <:AbstractVector{Vector{Int32}}
  cell_to_faces::Table{Int32,Vector{Int32},Vector{Int32}}
  cell_to_lface_to_pindex::Table{Int8,Vector{Int8},Vector{Int32}}
  cell_to_ctype::Vector{Int8}
  ctype_to_lface_to_own_ldofs::Vector{Vector{Vector{Int}}}
  ctype_to_num_dofs::Vector{Int}
  face_to_own_dofs::Table{Int32,Vector{Int32},Vector{Int32}}
  ctype_to_lface_to_pindex_to_pdofs::Vector{Vector{Vector{Vector{Int}}}}
end

Base.size(a::CellDofsNonOriented) = (length(a.cell_to_faces),)

Base.IndexStyle(::Type{<:CellDofsNonOriented}) = IndexLinear()

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
