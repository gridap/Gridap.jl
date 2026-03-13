
"""
    abstract type CellConformity <: GridapType

Minimum data required to describe dof ownership.
"""
abstract type CellConformity <: GridapType end

"""
"""
function get_cell_conformity(::SingleFieldFESpace)
  @abstractmethod
end

get_cell_conformity(f::UnconstrainedFESpace{V,<:CellConformity}) where V = f.metadata

"""
    struct GenericCellConformity <: CellConformity

A generic CellConformity implementation that stores cell-wise data.
"""
struct GenericCellConformity{A,B,C,D,E} <: CellConformity
  cell_lface_own_ldofs::A
  cell_lface_pindex_pdofs::B
  cell_d_num_dfaces::C
  cell_num_dofs::D
  cell_ldof_comp::E
end

Geometry.num_cells(c::GenericCellConformity) = length(c.cell_num_dofs)
Geometry.get_cell_type(c::GenericCellConformity) = Base.OneTo(num_cells(c))

"""
    GenericCellConformity(cell_lface_own_ldofs, cell_d_num_dfaces)
"""
function GenericCellConformity(
  cell_lface_own_ldofs :: AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}},
  cell_d_num_dfaces :: AbstractVector{<:AbstractVector{<:Integer}},
)
  cell_lface_pindex_pdofs = lazy_map(ReferenceFEs._trivial_face_own_dofs_permutations, cell_lface_own_ldofs)
  cell_num_dofs = lazy_map(x -> sum(length, x), cell_lface_own_ldofs)
  cell_ldof_comp = lazy_map(n -> fill(0,n), cell_num_dofs)
  GenericCellConformity(
    cell_lface_own_ldofs,
    cell_lface_pindex_pdofs,
    cell_d_num_dfaces,
    cell_num_dofs,
    cell_ldof_comp
  )
end

"""
    struct CompressedCellConformity{T} <: CellConformity

CellConformity implementation that compresses cell-wise data
on a few cell types.
"""
struct CompressedCellConformity{T} <: CellConformity
  cell_ctype::T
  ctype_lface_own_ldofs::Vector{Vector{Vector{Int}}}
  ctype_lface_pindex_pdofs::Vector{Vector{Vector{Vector{Int}}}}
  d_ctype_num_dfaces::Vector{Vector{Int}}
  ctype_num_dofs::Vector{Int}
  ctype_ldof_comp::Vector{Vector{Int}}
end

Geometry.num_cells(c::CompressedCellConformity) = length(c.cell_ctype)
Geometry.get_cell_type(c::CompressedCellConformity) = c.cell_ctype

function CellConformity(
  cell_ctype::AbstractVector{<:Integer},
  ctype_lface_own_ldofs,
  ctype_lface_pindex_pdofs,
  d_ctype_num_dfaces,
  ctype_num_dofs  = [sum(length, x) for x in ctype_lface_own_ldofs],
  ctype_ldof_comp = [fill(0,n) for n in ctype_num_dofs]
)
  CompressedCellConformity(
    cell_ctype,
    ctype_lface_own_ldofs,
    ctype_lface_pindex_pdofs,
    d_ctype_num_dfaces,
    ctype_num_dofs,
    ctype_ldof_comp
  )
end

"""
    CellConformity(cell_reffe::AbstractArray{<:ReferenceFE},conf::Conformity)
    CellConformity(cell_polys::AbstractArray{<:ReferenceFE},cell_basis::AbstractArray{<:AbstractArray{<:Field}},conf::Conformity)
"""
function CellConformity(cell_reffe::AbstractArray{<:ReferenceFE{D}},conf::Conformity) where D
  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
  ctype_lface_own_ldofs = map(reffe -> get_face_own_dofs(reffe,conf),ctype_reffe)
  ctype_lface_pindex_pdofs = map(reffe -> get_face_own_dofs_permutations(reffe,conf),ctype_reffe)
  d_ctype_num_dfaces = [ map(reffe -> num_faces(get_polytope(reffe),d),ctype_reffe) for d in 0:D]
  ctype_num_dofs = map(num_dofs, ctype_reffe)
  ctype_ldof_comp = map(get_dof_to_comp, ctype_reffe)
  return CellConformity(
    cell_ctype,ctype_lface_own_ldofs,ctype_lface_pindex_pdofs,d_ctype_num_dfaces,ctype_num_dofs,ctype_ldof_comp
  )
end

# function CellConformity(
#   cell_polys::AbstractArray{<:Polytope{D}},cell_basis::AbstractArray{<:AbstractArray{<:Field}},conf::Conformity
# ) where D
#   @check conf isa L2Conformity "Only L2 conformity is supported"
#   ctype_basis, cell_ctype = compress_cell_data(cell_basis)
#   ctype_poly, cell_ctype_poly = compress_cell_data(cell_polys)
#   @check isequal(cell_ctype,cell_ctype_poly) || isone(length(ctype_poly)) """
#     Inconsistent cell types between `cell_polys` and `cell_basis`
#   """
#   if isone(length(ctype_poly))
#     ctype_poly = fill(first(ctype_poly),length(ctype_basis))
#   end
# 
#   ctype_lface_own_ldofs = map(ctype_poly,ctype_basis) do p, b
#     nfaces = num_faces(p)
#     ndofs = length(b)
#     [ifelse(isequal(face,nfaces),collect(1:ndofs),Int[]) for face in 1:nfaces]
#   end
#   ctype_lface_pindex_pdofs = map(ReferenceFEs._trivial_face_own_dofs_permutations, ctype_lface_own_ldofs)
#   d_ctype_num_faces = [
#     map(p -> num_faces(p,d), ctype_poly) for d in 0:D
#   ]
#   ctype_num_dofs  = map(length, ctype_basis)
#   ctype_ldof_comp = [fill(0,n) for n in ctype_num_dofs]
# 
#   return CellConformity(
#     cell_ctype,ctype_lface_own_ldofs,ctype_lface_pindex_pdofs,d_ctype_num_faces,ctype_num_dofs,ctype_ldof_comp
#   )
# end

function Base.getproperty(a::CompressedCellConformity, sym::Symbol)
  if sym == :d_ctype_offset
    _d_ctype_offset(a)
  elseif sym == :d_ctype_ldface_own_ldofs
    _d_ctype_ldface_own_ldofs(a)
  else
    getfield(a, sym)
  end
end

function Base.propertynames(x::CompressedCellConformity, private::Bool=false)
  (fieldnames(typeof(x))...,:d_ctype_offset,:d_ctype_ldface_own_ldofs)
end

function _d_ctype_offset(a::CompressedCellConformity)
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

function _d_ctype_ldface_own_ldofs(a::CompressedCellConformity)
  num_ctypes = length(a.ctype_lface_own_ldofs)
  num_ds = length(a.d_ctype_num_dfaces)
  [[[ a.ctype_lface_own_ldofs[ctype][ldface+a.d_ctype_offset[d][ctype]]
    for ldface in 1:a.d_ctype_num_dfaces[d][ctype] ]
    for ctype in 1:num_ctypes ]
    for d in 1:num_ds ]
end

function get_d_ctype_lface_dofs(a::CompressedCellConformity, ctype_to_poly)
  ctype_lface_dofs = get_ctype_lface_dofs(a,ctype_to_poly)
  num_ctypes = length(a.ctype_lface_own_ldofs)
  num_ds = length(a.d_ctype_num_dfaces)
  [[[ ctype_lface_dofs[ctype][ldface+a.d_ctype_offset[d][ctype]]
    for ldface in 1:a.d_ctype_num_dfaces[d][ctype] ]
    for ctype in 1:num_ctypes ]
    for d in 1:num_ds ]
end

function get_ctype_lface_dofs(a::CompressedCellConformity, ctype_to_poly)
  return map(ReferenceFEs.face_own_data_to_face_data,ctype_to_poly,a.ctype_lface_own_ldofs)
end

"""
    struct CellFE <: GridapType

Minimum data required to build a conforming FE space.
At this moment, the some cell-wise info is compressed on cell types.
This can be relaxed in the future, and have an arbitrary cell-wise data.
"""
struct CellFE <: GridapType
  cell_conformity::CellConformity
  cell_shapefuns::AbstractArray{<:AbstractVector{<:Field}}
  cell_dof_basis::AbstractArray{<:AbstractVector{<:Dof}}
  cell_shapefuns_domain::DomainStyle
  cell_dof_basis_domain::DomainStyle
  max_order::Int
end
# If the shapefuns are not polynomials, max_order has to be understood as the order of a
# reasonable quadrature rule to integrate the shape functions. Only used by FESpace
# constructors that need to integrate the shape functions (e.g., ZeroMeanFESpace).

CellConformity(cell_fe::CellFE) = cell_fe.cell_conformity
Geometry.num_cells(cell_fe::CellFE) = num_cells(CellConformity(cell_fe))
Geometry.get_cell_type(cell_fe::CellFE) = get_cell_type(CellConformity(cell_fe))

"""
    CellFE(model::DiscreteModel,cell_reffe::AbstractArray{<:ReferenceFE},conformity::Conformity, args...)

Generate a CellFE from a vector of reference fes
"""
function CellFE(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:ReferenceFE},
  conformity::Conformity,
  args...
 )
  cell_conformity = CellConformity(cell_reffe,conformity)
  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
  cell_shapefuns = get_cell_shapefuns(model,cell_reffe,conformity,args...)
  cell_dof_basis = get_cell_dof_basis(model,cell_reffe,conformity,args...)
  cell_shapefuns_domain = ReferenceDomain()
  cell_dof_basis_domain = cell_shapefuns_domain
  max_order = maximum(map(get_order,ctype_reffe))
  CellFE(
    cell_conformity,
    cell_shapefuns,
    cell_dof_basis,
    cell_shapefuns_domain,
    cell_dof_basis_domain,
    max_order
  )
end

function get_cell_dof_basis(model::DiscreteModel,
                            cell_reffe::AbstractArray{<:ReferenceFE},
                            ::Conformity)
  lazy_map(get_dof_basis,cell_reffe)
end

function get_cell_shapefuns(model::DiscreteModel,
                            cell_reffe::AbstractArray{<:ReferenceFE},
                            ::Conformity)
  lazy_map(get_shapefuns,cell_reffe)
end

# Low level conforming FE Space constructor
# The user is not expected to call this function. Use the factory function instead
function _ConformingFESpace(
  vector_type::Type,
  model::DiscreteModel,
  face_labeling::FaceLabeling,
  cell_fe::CellFE,
  dirichlet_tags,
  dirichlet_components,
  trian = Triangulation(model))

  grid_topology = get_grid_topology(model)
  ntags = length(dirichlet_tags)

  cell_conformity = CellConformity(cell_fe)
  cell_dofs_ids, nfree, ndirichlet, dirichlet_dof_tag, dirichlet_cells = compute_conforming_cell_dofs(
    cell_conformity,grid_topology,face_labeling,dirichlet_tags,dirichlet_components
  )
  
  cell_shapefuns, cell_dof_basis = compute_cell_space(cell_fe,trian)

  cell_is_dirichlet = fill(false,num_cells(trian))
  cell_is_dirichlet[dirichlet_cells] .= true

  metadata = cell_conformity
  UnconstrainedFESpace(
    vector_type,
    nfree,
    ndirichlet,
    cell_dofs_ids,
    cell_shapefuns,
    cell_dof_basis,
    cell_is_dirichlet,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags,
    metadata
  )
end

function compute_cell_space(cell_shapefuns::AbstractArray{<:AbstractArray{<:Field}},
                            cell_dof_basis::AbstractArray{<:AbstractArray{<:Dof}},
                            cell_shapefuns_domain::DomainStyle,
                            cell_dof_basis_domain::DomainStyle,
                            trian::Triangulation)
  SingleFieldFEBasis(cell_shapefuns,trian,TestBasis(),cell_shapefuns_domain),
  CellDof(cell_dof_basis,trian,cell_dof_basis_domain)
end

function compute_cell_space(cellfe::CellFE,trian::Triangulation)
  compute_cell_space(
    cellfe.cell_shapefuns,
    cellfe.cell_dof_basis,
    cellfe.cell_shapefuns_domain,
    cellfe.cell_dof_basis_domain,
    trian
  )
end

function compute_cell_space(cell_fe,trian::Triangulation)
  cell_shapefuns, cell_dof_basis, d1, d2 = _compute_cell_space(cell_fe)
  SingleFieldFEBasis(cell_shapefuns,trian,TestBasis(),d1), CellDof(cell_dof_basis,trian,d2)
end

function _compute_cell_space(cell_fe::CellFE)
  ( cell_fe.cell_shapefuns,
    cell_fe.cell_dof_basis,
    cell_fe.cell_shapefuns_domain,
    cell_fe.cell_dof_basis_domain)
end

function _compute_cell_space(
  cell_reffe::AbstractArray{<:ReferenceFE})
  cell_shapefuns = lazy_map(get_shapefuns,cell_reffe)
  cell_dof_basis = lazy_map(get_dof_basis,cell_reffe)
  cell_shapefuns, cell_dof_basis, ReferenceDomain(), ReferenceDomain()
end

"""
The result is the tuple

    (cell_dofs, nfree, ndiri, dirichlet_dof_tag, dirichlet_cells)

If `dirichlet_components`  is given, then `get_dof_to_comp` has to be defined
for the reference elements in `reffes`.
"""
function compute_conforming_cell_dofs(
  cell_conformity::CompressedCellConformity,grid_topology,face_labeling,dirichlet_tags,dirichlet_components=nothing
)
  cell_to_ctype = cell_conformity.cell_ctype
  ctype_to_ldof_to_comp = cell_conformity.ctype_ldof_comp
  ctype_to_num_dofs = cell_conformity.ctype_num_dofs
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
    ctype_to_lface_to_pindex_to_pdofs) |> Table

  diri_cells = _generate_diri_cells(
    d_to_dface_to_tag,
    d_to_cell_to_dfaces)

  return cell_dofs, nfree, ndiri, diri_dof_tag, diri_cells
end

function _convert_dirichlet_components(dirichlet_tags::AbstractArray,dirichlet_components)
  dirichlet_components
end

function _convert_dirichlet_components(dirichlet_tags::Union{Integer,AbstractString},dirichlet_components)
  ifelse(isnothing(dirichlet_components),nothing,[dirichlet_components,])
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

"""
    generate_cell_dof_mask(
      scell_conformity::CellConformity,
      tcell_to_scell::AbstractVector{<:Integer},
      d_to_tcell_to_tdface,
      d_to_tdface_to_mask;
      reverse::Bool=false
    )

Given a `CellConformity` object, defined on a set of source cells (`scell`), and given 
a cell/face/edge/node masks on a set of target cells (`tcell`), this function generates a mask
for the degrees of freedom (dofs) on the target cells. 

Parameters:

- `scell_conformity`: The `CellConformity` object on the source cells (`scell`).
- `tcell_to_scell`: A vector mapping target cells (`tcell`) to source cells (`scell`).
- `d_to_tcell_to_tdface`: For each dimension `d`, a `Table` mapping target cells to its d-faces (edges, faces, ...)
- `d_to_tdface_to_mask`: For each dimension `d`, an array of booleans indicating whether the d-face is masked or not.

Returns:

- `tcell_dof_mask`: A vector that for each target cell contains a boolean vector of size equal
   to the number of dofs in that cell, containing the mask.

Modes of operation:

- If `reverse = false`, the function generates a mask where the dofs are
  `true` if the corresponding d-face is masked (`true`).
- If `reverse = true`, the mask is reversed, meaning that the dofs are `true` if 
  the corresponding d-face is not masked (`false`).
"""
function generate_cell_dof_mask(
  scell_conformity::CompressedCellConformity,
  tcell_to_scell::AbstractVector{<:Integer},
  d_to_tcell_to_tdface,
  d_to_tdface_to_mask;
  reverse::Bool=false
)
  scell_ctype = scell_conformity.cell_ctype
  ctype_ndofs = scell_conformity.ctype_num_dofs
  d_ctype_ldface_own_ldofs = scell_conformity.d_ctype_ldface_own_ldofs
  D = length(d_ctype_ldface_own_ldofs) - 1

  caches = map(array_cache, d_to_tcell_to_tdface)
  tcell_dof_mask = Vector{Vector{Bool}}(undef, length(tcell_to_scell))
  for (tcell,scell) in enumerate(tcell_to_scell)
    if scell < 1
      tcell_dof_mask[tcell] = Bool[]
      continue
    end

    ctype = scell_ctype[scell]
    dof_mask = fill(false,ctype_ndofs[ctype])

    for d in 0:D
      ldface_to_ldofs = d_ctype_ldface_own_ldofs[d+1][ctype]
      tdfaces = getindex!(caches[d+1],d_to_tcell_to_tdface[d+1],tcell)
      for (ldface,tdface) in enumerate(tdfaces)
        if d_to_tdface_to_mask[d+1][tdface]
          ldofs = ldface_to_ldofs[ldface]
          dof_mask[ldofs] .= true
        end
      end
    end

    tcell_dof_mask[tcell] = xor.(dof_mask,reverse)
  end

  return tcell_dof_mask
end

function generate_dof_mask(
  scell_conformity::CompressedCellConformity,
  scell_dof_ids::AbstractVector,
  tcell_to_scell::AbstractVector{<:Integer},
  d_to_tcell_to_tdface,
  d_to_tdface_to_mask,
  n_free_dofs;
  reverse::Bool=false
)
  scell_ctype = get_cell_type(scell_conformity)
  ctype_lface_own_ldofs = scell_conformity.ctype_lface_own_ldofs
  d_ctype_offset = scell_conformity.d_ctype_offset
  D = length(d_ctype_offset) - 1

  caches = map(array_cache, d_to_tcell_to_tdface)
  dofs_cache = array_cache(scell_dof_ids)

  dof_to_mask = fill(false,n_free_dofs)
  for (tcell,scell) in enumerate(tcell_to_scell)
    (scell < 1) && continue

    ctype = scell_ctype[scell]
    dofs = getindex!(dofs_cache,scell_dof_ids,scell)
    lface_own_ldofs = ctype_lface_own_ldofs[ctype]

    for d in 0:D
      offset = d_ctype_offset[d+1][ctype]
      dfaces = getindex!(caches[d+1],d_to_tcell_to_tdface[d+1],tcell)
      for (ldface,dface) in enumerate(dfaces)
        mask = d_to_tdface_to_mask[d+1][dface]
        ldofs = lface_own_ldofs[offset+ldface]
        for ldof in ldofs
          dof = dofs[ldof]
          if dof > 0 # Avoid dirichlet dofs
            dof_to_mask[dof] |= mask
          end
        end
      end
    end
  end

  if reverse
    dof_to_mask .= .!dof_to_mask
  end
  return dof_to_mask
end

"""
    generate_dof_mask(
      space::FESpace, labels::FaceLabeling, tags; 
      cell_conformity = get_cell_conformity(space), reverse::Bool=false
    )

Generate a mask for the dofs in the FESpace `space`, based on the face labeling `labels`
and the `tags` provided. If `reverse` is `false`, the dofs are masked (`true`) if they 
belong to a tagged face. If `reverse` is `true`, the mask is reversed.
"""
function generate_dof_mask(
  space::FESpace, labels::FaceLabeling, tags; 
  cell_conformity = get_cell_conformity(space), reverse::Bool=false
)
  trian = get_triangulation(space)
  model = get_background_model(trian)
  topo  = get_grid_topology(model)

  Df = num_cell_dims(trian)
  mface_to_tface = get_glue(trian,Val(Df)).mface_to_tface
  cell_dof_ids = get_cell_dof_ids(space)
  d_to_dface_to_mask = [get_face_mask(labels,tags,d) for d in 0:Df]
  d_to_face_to_dface = [Geometry.get_faces(topo,Df,d) for d in 0:Df]

  return generate_dof_mask(
    cell_conformity, cell_dof_ids, mface_to_tface,
    d_to_face_to_dface, d_to_dface_to_mask, num_free_dofs(space);
    reverse
  )
end

struct CellDofsNonOriented <:AbstractVector{Vector{Int32}}
  cell_to_faces::Table{Int32,Vector{Int32},Vector{Int32}}
  cell_to_lface_to_pindex::Table{GridapLocalInt,Vector{GridapLocalInt},Vector{Int32}}
  cell_to_ctype::Vector{Int32}
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
