
"""
    generate_cells_around(cell_to_faces::Table)
    generate_cells_around(cell_to_faces::Table, nfaces::Int)
"""
function generate_cells_around(
  cell_to_faces::Table,
  nfaces::Int = maximum(cell_to_faces.data))

  data, ptrs = _face_to_cells(cell_to_faces.data,cell_to_faces.ptrs,nfaces)
  Table(data,ptrs)
end

"""
    generate_cell_to_faces(
      cell_to_vertices::Table,
      cell_type_to_lface_to_lvertices::Vector{Vector{Vector{Int}}},
      cell_to_cell_type::AbstractVector{<:Integer},
      vertex_to_cells::Table) -> Table
"""
function generate_cell_to_faces(
  cell_to_vertices::Table,
  cell_type_to_lface_to_lvertices::Vector{Vector{Vector{Int}}},
  cell_to_cell_type::AbstractVector{<:Integer},
  vertex_to_cells::Table)

  data, ptrs = _generate_cell_to_faces(
    cell_to_vertices.data,
    cell_to_vertices.ptrs,
    cell_type_to_lface_to_lvertices,
    cell_to_cell_type,
    vertex_to_cells.data,
    vertex_to_cells.ptrs)

  Table(data,ptrs)
end

"""
    generate_face_to_face_type(
      cell_to_faces::Table,
      cell_to_cell_type::AbstractVector{<:Integer},
      cell_type_to_lface_to_face_type::Vector{Vector{T}} [,
      nfaces::Int ]) where T<:Integer -> Vector{T}
"""
function generate_face_to_face_type(
  cell_to_faces::Table,
  cell_to_cell_type::AbstractVector{<:Integer},
  cell_type_to_lface_to_face_type::Vector{Vector{T}},
  nfaces::Int=maximum(cell_to_faces.data)) where T<:Integer

  _generate_face_to_ftype(
    cell_to_faces_data,
    cell_to_faces_ptrs,
    cell_to_cell_type,
    cell_type_to_lface_to_face_type,
    nfaces)

end



function find_cell_to_faces(
  cell_to_vertices_data::AbstractVector{<:Integer},
  cell_to_vertices_ptrs::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices,
  cell_to_ctype,
  vertex_to_faces_data::AbstractVector{<:Integer},
  vertex_to_faces_ptrs::AbstractVector{<:Integer})
  _cell_to_faces_from_vertex_to_faces(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices,
    cell_to_ctype,
    vertex_to_faces_data,
    vertex_to_faces_ptrs)
end

function generate_facet_to_isboundary(face_to_cells_ptrs::AbstractVector{<:Integer})
  nfaces = length(face_to_cells_ptrs)-1
  face_to_isboundary = fill(false,nfaces)
  _generate_face_to_isboundary_fill!(face_to_isboundary,face_to_cells_ptrs)
  face_to_isboundary
end

function generate_face_to_isboundary(
  facet_to_isboundary::AbstractVector{Bool},
  face_to_facets_data::AbstractVector{<:Integer},
  face_to_facets_ptrs::AbstractVector{<:Integer})
  nobjects = length(face_to_facets_ptrs)-1
  object_to_isboundary = fill(false,nobjects)
  _generate_object_to_isboundary_fill!(
    object_to_isboundary,
    facet_to_isboundary,
    face_to_facets_data,
    face_to_facets_ptrs)
  object_to_isboundary
end

function generate_face_to_vertices(
  cell_to_vertices_data::AbstractVector{<:Integer},
  cell_to_vertices_ptrs::AbstractVector{<:Integer},
  cell_to_faces_data::AbstractVector{<:Integer},
  cell_to_faces_ptrs::AbstractVector{<:Integer},
  cell_to_ctype::AbstractVector{<:Integer},
  ctype_to_lface_to_lvertices::Vector{Vector{Vector{Int}}},
  nfaces::Int=maximum(cell_to_faces_data))
  _face_to_vertices(
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    cell_to_faces_data,
    cell_to_faces_ptrs,
    cell_to_ctype,
    ctype_to_lface_to_lvertices,
    nfaces)
end

function refine_grid_connectivity(
  cell_to_points_data::AbstractVector{T},
  cell_to_points_ptrs::AbstractVector{P},
  ltcell_to_lpoints) where {T,P}

  nltcells = length(ltcell_to_lpoints)
  ncells = length(cell_to_points_ptrs) - 1
  ntcells = ncells * nltcells

  tcell_to_points_ptrs = zeros(P,ntcells+1)

  _refine_grid_connectivity_count!(
    tcell_to_points_ptrs,
    ncells,
    ltcell_to_lpoints)

  length_to_ptrs!(tcell_to_points_ptrs)

  ndata = tcell_to_points_ptrs[end]-1

  tcell_to_points_data = zeros(T,ndata)

  _refine_grid_connectivity!(
    tcell_to_points_data,
    cell_to_points_data,
    cell_to_points_ptrs,
    ltcell_to_lpoints )

  (tcell_to_points_data, tcell_to_points_ptrs)

end

function generate_tface_to_face(
  cell_to_faces_data::AbstractVector{T},
  cell_to_faces_ptrs,
  tcell_to_tfaces_data,
  tcell_to_tfaces_ptrs,
  ltcell_to_lnodes,
  ltface_to_ltnodes,
  lface_to_lnodes,
  ntfaces) where T

  ltcell_to_lfaces =  _generate_local_face_relations(
    ltcell_to_lnodes,
    ltface_to_ltnodes,
    lface_to_lnodes)

  _generate_tface_to_face(
    cell_to_faces_data,
    cell_to_faces_ptrs,
    tcell_to_tfaces_data,
    tcell_to_tfaces_ptrs,
    ltcell_to_lfaces,
    ntfaces)

end

function find_gface_to_face(
  face_to_nodes_data,
  face_to_nodes_ptrs,
  node_to_faces_data::AbstractVector{T},
  node_to_faces_ptrs,
  gface_to_nodes_data,
  gface_to_nodes_ptrs) where T

  ngfaces = length(gface_to_nodes_ptrs) - 1
  gface_to_face = zeros(T,ngfaces)
  n = max_cells_arround_vertex(node_to_faces_ptrs)
  faces_around = fill(UNSET,n)
  faces_around_scratch = fill(UNSET,n)

  _fill_gface_to_face!(
    gface_to_face,
    face_to_nodes_data,
    face_to_nodes_ptrs,
    node_to_faces_data,
    node_to_faces_ptrs,
    gface_to_nodes_data,
    gface_to_nodes_ptrs,
    faces_around,
    faces_around_scratch)

  gface_to_face

end

# Helpers

const UNSET = 0

function _generate_face_to_isboundary_fill!(
  face_to_isboundary, face_to_cells_ptrs)
  nfaces = length(face_to_isboundary)
  for face in 1:nfaces
    ncells_around = face_to_cells_ptrs[face+1] - face_to_cells_ptrs[face]
    if ncells_around == 1
      face_to_isboundary[face] = true
    elseif ncells_around == 2
      face_to_isboundary[face] = false
    else
      @unreachable
    end
  end
end

function  _generate_object_to_isboundary_fill!(
  object_to_isboundary,
  face_to_isboundary,
  object_to_faces_data,
  object_to_faces_ptrs)

  nobjects = length(object_to_isboundary)
  for object in 1:nobjects
    a = object_to_faces_ptrs[object]-1
    b = object_to_faces_ptrs[object+1]
    nlfaces = b-(a+1)
    for lface in 1:nlfaces
      face = object_to_faces_data[a+lface]
      isboundary = face_to_isboundary[face]
      if isboundary
        object_to_isboundary[object] = true
        break
      end
    end
  end

end

function _face_to_cells(
  cell_to_faces_data::AbstractVector{L},
  cell_to_faces_ptrs::AbstractVector{P},
  nfaces) where {L,P}

  face_to_cells_ptrs = zeros(P,nfaces+1)

  _face_to_cells_count!(
    face_to_cells_ptrs, cell_to_faces_data, cell_to_faces_ptrs)

  length_to_ptrs!(face_to_cells_ptrs)

  ndata = face_to_cells_ptrs[end]-1

  face_to_cells_data = Vector{L}(undef,ndata)

  _face_to_cells_fill!(
    face_to_cells_data, face_to_cells_ptrs,
    cell_to_faces_data, cell_to_faces_ptrs)

  rewind_ptrs!(face_to_cells_ptrs)

  (face_to_cells_data, face_to_cells_ptrs)

end

function _face_to_cells_count!(
    face_to_cells_ptrs, cell_to_faces_data, cell_to_faces_ptrs)

  ncells = length(cell_to_faces_ptrs) - 1
  for cell in 1:ncells
    a = cell_to_faces_ptrs[cell]
    b = cell_to_faces_ptrs[cell+1]-1
    for p in a:b
      face = cell_to_faces_data[p]
      face_to_cells_ptrs[face+1] += 1
    end
  end

end

function _face_to_cells_fill!(
    face_to_cells_data, face_to_cells_ptrs,
    cell_to_faces_data, cell_to_faces_ptrs)

  ncells = length(cell_to_faces_ptrs) - 1
  for cell in 1:ncells
    a = cell_to_faces_ptrs[cell]
    b = cell_to_faces_ptrs[cell+1]-1
    for p in a:b
      face = cell_to_faces_data[p]
      q = face_to_cells_ptrs[face]
      face_to_cells_data[q] = cell
      face_to_cells_ptrs[face] += 1
    end
  end

end

function _generate_cell_to_faces(
    cell_to_vertices_data::AbstractVector{L},
    cell_to_vertices_ptrs::AbstractVector{P},
    ctype_to_lface_to_lvertices,
    cell_to_ctype,
    vertex_to_cells_data,
    vertex_to_cells_ptrs) where {L,P}

  cell_to_faces_ptrs = _cell_to_faces_count(
    cell_to_vertices_ptrs,
    cell_to_ctype,
    ctype_to_lface_to_lvertices)

  length_to_ptrs!(cell_to_faces_ptrs)

  ndata = cell_to_faces_ptrs[end]-1

  cell_to_faces_data = fill(L(UNSET),ndata)

  nvertices = max_nvertices_in_lface(ctype_to_lface_to_lvertices)
  vertices = fill(UNSET,nvertices)
  vertices_scratch = fill(UNSET,nvertices)

  nvertices = max_cells_arround_vertex(vertex_to_cells_ptrs)
  cells_around = fill(UNSET,nvertices)
  cells_around_scratch = fill(UNSET,nvertices)
  lface_of_cells_around = fill(UNSET,nvertices)

  _cell_to_faces_fill!(
    cell_to_faces_data,
    cell_to_faces_ptrs,
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices,
    cell_to_ctype,
    vertex_to_cells_data,
    vertex_to_cells_ptrs,
    vertices,
    vertices_scratch,
    cells_around,
    cells_around_scratch,
    lface_of_cells_around)

  (cell_to_faces_data, cell_to_faces_ptrs)

end

function  _cell_to_faces_from_vertex_to_faces(
  cell_to_vertices_data::AbstractVector{L},
  cell_to_vertices_ptrs::AbstractVector{P},
  ctype_to_lface_to_lvertices,
  cell_to_ctype,
  vertex_to_faces_data,
  vertex_to_faces_ptrs) where {L,P}

  cell_to_faces_ptrs = _cell_to_faces_count(
    cell_to_vertices_ptrs,
    cell_to_ctype,
    ctype_to_lface_to_lvertices)

  length_to_ptrs!(cell_to_faces_ptrs)

  ndata = cell_to_faces_ptrs[end]-1

  cell_to_faces_data = fill(L(UNSET),ndata)

  nvertices = max_nvertices_in_lface(ctype_to_lface_to_lvertices)
  vertices = fill(UNSET,nvertices)

  nvertices = max_cells_arround_vertex(vertex_to_faces_ptrs)
  faces_around = fill(UNSET,nvertices)
  faces_around_scratch = fill(UNSET,nvertices)

  _cell_to_faces_fill!(
    cell_to_faces_ptrs,
    cell_to_faces_data,
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices,
    cell_to_ctype,
    vertex_to_faces_data,
    vertex_to_faces_ptrs,
    vertices,
    faces_around,
    faces_around_scratch)

  (cell_to_faces_data, cell_to_faces_ptrs)
end

function _cell_to_faces_count(
  cell_to_vertices_ptrs::AbstractVector{P},
  cell_to_ctype,
  ctype_to_lface_to_lvertices) where P

  ncells = length(cell_to_vertices_ptrs) - 1

  cell_to_faces_ptrs = zeros(P,ncells+1)

  type_to_nlfaces = _type_to_nlfaces(ctype_to_lface_to_lvertices)

  _cell_to_faces_count!(
    cell_to_faces_ptrs,
    type_to_nlfaces,
    cell_to_ctype)

  cell_to_faces_ptrs

end

function _cell_to_faces_fill!(
  cell_to_faces_ptrs,
  cell_to_faces_data,
  cell_to_vertices_data,
  cell_to_vertices_ptrs,
  ctype_to_lface_to_lvertices,
  cell_to_ctype,
  vertex_to_faces_data,
  vertex_to_faces_ptrs,
  vertices,
  faces_around,
  faces_around_scratch)

  ncells = length(cell_to_ctype)

  for cell in 1:ncells

    ctype = cell_to_ctype[cell]
    lface_to_lvertices = ctype_to_lface_to_lvertices[ctype]
    nlfaces = length(lface_to_lvertices)
    a = cell_to_faces_ptrs[cell]-1

    for lface in 1:nlfaces

      _fill_vertices_in_lface!(
        vertices,
        lface,
        lface_to_lvertices,
        cell,
        cell_to_vertices_data,
        cell_to_vertices_ptrs)

      _find_cells_around_vertices!(
        faces_around,
        faces_around_scratch,
        vertices,
        vertex_to_faces_data,
        vertex_to_faces_ptrs)

      for face in faces_around
        if face != UNSET
          cell_to_faces_data[a+lface] = face
          break
        end
      end

    end

  end

end

function max_nvertices_in_lface(ctype_to_lface_to_lvertices)
  n = 0
  for lface_to_lvertices in ctype_to_lface_to_lvertices
    for lvertices in lface_to_lvertices
      m = length(lvertices)
      n = max(n,m)
    end
  end
  n
end

function max_cells_arround_vertex(vertex_to_cells_ptrs)
  n = 0
  nvertices = length(vertex_to_cells_ptrs)-1
  for vertex in 1:nvertices
    m = vertex_to_cells_ptrs[vertex+1]-vertex_to_cells_ptrs[vertex]
    n = max(n,m)
  end
  n
end

function _type_to_nlfaces(ctype_to_lface_to_lvertices)
  [ length(lface_to_lvertices)
   for lface_to_lvertices in ctype_to_lface_to_lvertices]
end

function _cell_to_faces_count!(
    cell_to_faces_ptrs,
    type_to_nlfaces,
    cell_to_ctype)

  ncells = length(cell_to_ctype)
  for cell in 1:ncells
    ctype = cell_to_ctype[cell]
    nlfaces = type_to_nlfaces[ctype]
    cell_to_faces_ptrs[cell+1] = nlfaces
  end

end

function  _cell_to_faces_fill!(
    cell_to_faces_data,
    cell_to_faces_ptrs,
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    ctype_to_lface_to_lvertices,
    cell_to_ctype,
    vertex_to_cells_data,
    vertex_to_cells_ptrs,
    vertices,
    vertices_scratch,
    cells_around,
    cells_around_scratch,
    lface_of_cells_around)

  face = 1

  ncells = length(cell_to_ctype)

  for cell in 1:ncells

    ctype = cell_to_ctype[cell]
    lface_to_lvertices = ctype_to_lface_to_lvertices[ctype]
    nlfaces = length(lface_to_lvertices)
    a = cell_to_faces_ptrs[cell]-1

    for lface in 1:nlfaces

      if cell_to_faces_data[a+lface] != UNSET
        continue
      end

      _fill_vertices_in_lface!(
        vertices,
        lface,
        lface_to_lvertices,
        cell,
        cell_to_vertices_data,
        cell_to_vertices_ptrs)

      _find_cells_around_vertices!(
        cells_around,
        cells_around_scratch,
        vertices,
        vertex_to_cells_data,
        vertex_to_cells_ptrs)

      _find_lface_of_cells_around_lface!(
        lface_of_cells_around,
        cells_around,
        lface,
        ctype_to_lface_to_lvertices,
        cell,
        cell_to_vertices_data,
        cell_to_vertices_ptrs,
        cell_to_ctype,
        vertices,
        vertices_scratch)

      _fill_face_in_cells_arround!(
        cell_to_faces_data,
        cell_to_faces_ptrs,
        face,
        cells_around,
        lface_of_cells_around)

      face += 1
    end

  end

end

function _fill_face_in_cells_arround!(
  cell_to_faces_data,
  cell_to_faces_ptrs,
  face,
  cells_around,
  lface_of_cells_around)

  for icell_around in 1:length(cells_around)
    cell_around = cells_around[icell_around]
    if cell_around == UNSET
      continue
    end
    lface = lface_of_cells_around[icell_around]
    f = cell_to_faces_ptrs[cell_around]-1
    cell_to_faces_data[f+lface] = face
  end

end

function _find_cells_around_vertices!(
  cells_around,
  cells_around_scratch,
  vertices,
  vertex_to_cells_data,
  vertex_to_cells_ptrs)

  ncells_around = UNSET
  ncells_around_scratch = UNSET

  for ivertex in 1:length(vertices)
    vertex = vertices[ivertex]
    if vertex == UNSET
      continue
    end
    if ivertex == 1
      ncells_around = _fill_cells_around_scratch!(
        cells_around,
        vertex,
        vertex_to_cells_data,
        vertex_to_cells_ptrs)
    else
      ncells_around_scratch = _fill_cells_around_scratch!(
        cells_around_scratch,
        vertex,
        vertex_to_cells_data,
        vertex_to_cells_ptrs)
      _set_intersection!(
        cells_around,cells_around_scratch,
        ncells_around,ncells_around_scratch)
    end
  end

end

function _fill_cells_around_scratch!(
  cells_around_scratch,
  vertex,
  vertex_to_cells_data,
  vertex_to_cells_ptrs)

  cells_around_scratch .= UNSET
  d = vertex_to_cells_ptrs[vertex] - 1
  ncells_around = vertex_to_cells_ptrs[vertex+1] - (d + 1)
  for icell_around in 1:ncells_around
    cell_around = vertex_to_cells_data[d+icell_around]
    cells_around_scratch[icell_around] = cell_around
  end
  ncells_around

end

function _set_intersection!(
  cells_around,cells_around_scratch,
  ncells_around,ncells_around_scratch)
  for i in 1:ncells_around
    if cells_around[i] == UNSET
      continue
    end
    _find_eq!(i,cells_around,cells_around_scratch,ncells_around_scratch)
  end
end

function _find_eq!(i,cells_around,cells_around_scratch,ncells_around_scratch)
  for j in 1:ncells_around_scratch
    if cells_around[i] == cells_around_scratch[j]
      return
    end
  end
  cells_around[i] = UNSET
  return
end

function _find_lface_of_cells_around_lface!(
  lface_of_cells_around,
  cells_around,
  lface,
  ctype_to_lface_to_lvertices,
  cell,
  cell_to_vertices_data,
  cell_to_vertices_ptrs,
  cell_to_ctype,
  vertices,
  vertices_scratch)

  lface_of_cells_around .= UNSET

  for icell_around in 1:length(cells_around)
    cell_around = cells_around[icell_around]
    if cell_around == UNSET
      continue
    end

    lface_of_cell_around = _find_lface_with_same_vertices(
      vertices,
      vertices_scratch,
      cell_around,
      cell_to_vertices_data,
      cell_to_vertices_ptrs,
      cell_to_ctype,
      ctype_to_lface_to_lvertices)

    lface_of_cells_around[icell_around] = lface_of_cell_around

  end

end

function _find_lface_with_same_vertices(
  vertices,
  vertices_scratch,
  cell,
  cell_to_vertices_data,
  cell_to_vertices_ptrs,
  cell_to_ctype,
  ctype_to_lface_to_lvertices)

  ctype = cell_to_ctype[cell]
  lface_to_lvertices = ctype_to_lface_to_lvertices[ctype]

  nlfaces = length(lface_to_lvertices)
  c = cell_to_vertices_ptrs[cell]-1

  for lface in 1:nlfaces

    _fill_vertices_in_lface!(
      vertices_scratch,
      lface,
      lface_to_lvertices,
      cell,
      cell_to_vertices_data,
      cell_to_vertices_ptrs)

    if _set_equal(vertices,vertices_scratch)
      return lface
    end

  end

  return UNSET

end

function _fill_vertices_in_lface!(
  vertices,
  lface,
  lface_to_lvertices,
  cell,
  cell_to_vertices_data,
  cell_to_vertices_ptrs)

  vertices .= UNSET
  c = cell_to_vertices_ptrs[cell] - 1
  lvertices = lface_to_lvertices[lface]
  nlfvertex = length(lvertices)
  for lfvertex in 1:nlfvertex
    lvertex = lvertices[lfvertex]
    vertex = cell_to_vertices_data[c+lvertex]
    vertices[lfvertex] = vertex
  end

end

function _set_equal(vertices,vertices_scratch)

  b = _is_subset(vertices,vertices_scratch)
  if b == false; return false; end
  b = _is_subset(vertices_scratch,vertices)
  if b == false; return false; end
  return true

end

function _is_subset(vertices,vertices_scratch)
  for i in 1:length(vertices)
    v = vertices[i]
    if v == UNSET
      continue
    end
    b = _find_eq(v,vertices_scratch)
    if b == false; return false; end
  end
  return true
end

function _find_eq(v,vertices_scratch)
  for vs in vertices_scratch
    if v == vs
      return true
    end
  end
  return false
end

function _generate_face_to_ftype(
  cell_to_faces_data,
  cell_to_faces_ptrs,
  cell_to_ctype,
  ctype_to_lface_to_ftype::Vector{Vector{T}},
  nfaces) where T

  face_to_ftype = fill(T(UNSET),nfaces)
  _face_to_ftype_fill!(
    face_to_ftype,
    cell_to_faces_data,
    cell_to_faces_ptrs,
    cell_to_ctype,
    ctype_to_lface_to_ftype)
  face_to_ftype
end

function _face_to_ftype_fill!(
  face_to_ftype,
  cell_to_faces_data,
  cell_to_faces_ptrs,
  cell_to_ctype,
  ctype_to_lface_to_ftype)

  ncells = length(cell_to_ctype)
  for cell in 1:ncells
    ctype = cell_to_ctype[cell]
    lface_to_ftype = ctype_to_lface_to_ftype[ctype]
    a = cell_to_faces_ptrs[cell]-1
    nlfaces = cell_to_faces_ptrs[cell+1] - (a + 1)
    for lface in 1:nlfaces
      face = cell_to_faces_data[a+lface]
      if face_to_ftype[face] != UNSET
        continue
      end
      ftype = lface_to_ftype[lface]
      face_to_ftype[face] = ftype
    end
  end

end

function _face_to_vertices(
  cell_to_vertices_data::AbstractVector{L},
  cell_to_vertices_ptrs::AbstractVector{P},
  cell_to_faces_data,
  cell_to_faces_ptrs,
  cell_to_ctype,
  ctype_to_lface_to_lvertices,
  nfaces) where {L,P}

  face_to_vertices_ptrs = fill(P(UNSET),nfaces+1)

  _face_to_vertices_count!(
    face_to_vertices_ptrs,
    cell_to_faces_data,
    cell_to_faces_ptrs,
    cell_to_ctype,
    ctype_to_lface_to_lvertices)

  length_to_ptrs!(face_to_vertices_ptrs)
  ndata = face_to_vertices_ptrs[end]-1
  face_to_vertices_data = fill(L(UNSET),ndata)

  _face_to_vertices_fill!(
    face_to_vertices_data,
    face_to_vertices_ptrs,
    cell_to_vertices_data,
    cell_to_vertices_ptrs,
    cell_to_faces_data,
    cell_to_faces_ptrs,
    cell_to_ctype,
    ctype_to_lface_to_lvertices)

  (face_to_vertices_data, face_to_vertices_ptrs)

end

function _face_to_vertices_count!(
  face_to_vertices_ptrs,
  cell_to_faces_data,
  cell_to_faces_ptrs,
  cell_to_ctype,
  ctype_to_lface_to_lvertices)

  ncells = length(cell_to_ctype)
  for cell in 1:ncells
    ctype = cell_to_ctype[cell]
    lface_to_lvertices = ctype_to_lface_to_lvertices[ctype]
    a = cell_to_faces_ptrs[cell]-1
    nlfaces = cell_to_faces_ptrs[cell+1] - (a + 1)
    for lface in 1:nlfaces
      face = cell_to_faces_data[a+lface]
      if face_to_vertices_ptrs[face+1] != UNSET
        continue
      end
      nfvertices = length(lface_to_lvertices[lface])
      face_to_vertices_ptrs[face+1] = nfvertices
    end
  end

end

function _face_to_vertices_fill!(
  face_to_vertices_data,
  face_to_vertices_ptrs,
  cell_to_vertices_data,
  cell_to_vertices_ptrs,
  cell_to_faces_data,
  cell_to_faces_ptrs,
  cell_to_ctype,
  ctype_to_lface_to_lvertices)

  ncells = length(cell_to_ctype)
  for cell in 1:ncells
    ctype = cell_to_ctype[cell]
    lface_to_lvertices = ctype_to_lface_to_lvertices[ctype]
    a = cell_to_faces_ptrs[cell]-1
    c = cell_to_vertices_ptrs[cell]-1
    nlfaces = cell_to_faces_ptrs[cell+1] - (a + 1)
    for lface in 1:nlfaces
      face = cell_to_faces_data[a+lface]
      v = face_to_vertices_ptrs[face]-1
      if face_to_vertices_data[v+1] != UNSET
        continue
      end
      lvertices = lface_to_lvertices[lface]
      nfvertices = length(lvertices)
      for lfvertex in 1:nfvertices
        lvertex = lvertices[lfvertex]
        vertex = cell_to_vertices_data[c+lvertex]
        face_to_vertices_data[v+lfvertex] = vertex
      end
    end
  end

end

function  _refine_grid_connectivity_count!(
    tcell_to_points_ptrs,
    ncells,
    ltcell_to_lpoints)

  tcell = 1

  for cell in 1:ncells
    for lpoints in ltcell_to_lpoints
      tcell_to_points_ptrs[tcell+1] = length(lpoints)
      tcell +=1
    end
  end

end

function _refine_grid_connectivity!(
    tcell_to_points_data,
    cell_to_points_data,
    cell_to_points_ptrs,
    ltcell_to_lpoints )

  ncells = length(cell_to_points_ptrs) - 1

  k = 1
  for cell in 1:ncells
    a = cell_to_points_ptrs[cell]-1
    for lpoints in ltcell_to_lpoints
      for lpoint in lpoints
        point = cell_to_points_data[a+lpoint]
        tcell_to_points_data[k] = point
        k += 1
      end
    end
  end

end

function _generate_local_face_relations(
  ltcell_to_lnodes, ltface_to_ltnodes, lface_to_lnodes)

  nltcells = length(ltcell_to_lnodes)
  nltfaces = length(ltface_to_ltnodes)
  nlfaces = length(lface_to_lnodes)

  ltcell_ltface_to_lface = [ zeros(Int,nltfaces) for ltcell in 1:nltcells  ]

  for ltcell in 1:nltcells

    ltnode_to_lnode = ltcell_to_lnodes[ltcell]

    for ltface in 1:nltfaces
      ltnodes = ltface_to_ltnodes[ltface]
      for lface in 1:nlfaces
        lnodes = lface_to_lnodes[lface]
        allin = true
        for ltnode in ltnodes
          lnode = ltnode_to_lnode[ltnode]
          if !(lnode in lnodes)
            allin = false
            break
          end
        end
        if allin
          ltcell_ltface_to_lface[ltcell][ltface] = lface
          break
        end
      end
    end

  end

  ltcell_ltface_to_lface

end


function _generate_tface_to_face(
  cell_to_faces_data::AbstractVector{T},
  cell_to_faces_ptrs,
  tcell_to_tfaces_data,
  tcell_to_tfaces_ptrs,
  ltcell_to_lfaces,
  ntfaces) where T

  tface_to_face = zeros(T,ntfaces)

  _generate_tface_to_face!(
    tface_to_face,
    cell_to_faces_data,
    cell_to_faces_ptrs,
    tcell_to_tfaces_data,
    tcell_to_tfaces_ptrs,
    ltcell_to_lfaces)

  tface_to_face

end

function  _generate_tface_to_face!(
  tface_to_face,
  cell_to_faces_data,
  cell_to_faces_ptrs,
  tcell_to_tfaces_data,
  tcell_to_tfaces_ptrs,
  ltcell_to_lfaces)

  ncells = length(cell_to_faces_ptrs)-1
  nltcells = length(ltcell_to_lfaces)

  tcell = 1
  for cell in 1:ncells
    c = cell_to_faces_ptrs[cell]-1
    for ltcell in 1:nltcells
      a = tcell_to_tfaces_ptrs[tcell]-1
      b = tcell_to_tfaces_ptrs[tcell+1]
      nltfaces = b - (a+1)
      ltface_to_lface = ltcell_to_lfaces[ltcell]
      for ltface in 1:nltfaces
        tface = tcell_to_tfaces_data[a+ltface]
        lface = ltface_to_lface[ltface]
        if lface != 0
          face = cell_to_faces_data[c+lface]
          tface_to_face[tface] = face
        end
      end
      tcell += 1
    end
  end

end

function  _fill_gface_to_face!(
  gface_to_face,
  face_to_nodes_data,
  face_to_nodes_ptrs,
  node_to_faces_data,
  node_to_faces_ptrs,
  gface_to_nodes_data,
  gface_to_nodes_ptrs,
  faces_around,
  faces_around_scratch)

  ngfaces = length(gface_to_nodes_ptrs) - 1

  nfaces_around = UNSET
  nfaces_around_scratch = UNSET

  for gface in 1:ngfaces

    a = gface_to_nodes_ptrs[gface]-1
    b = gface_to_nodes_ptrs[gface+1]
    nlnodes = b-(a+1)

    for lnode in 1:nlnodes
      node = gface_to_nodes_data[lnode+a]
      if lnode == 1
        nfaces_around = _fill_cells_around_scratch!(
          faces_around,
          node,
          node_to_faces_data,
          node_to_faces_ptrs)
      else
        nfaces_around_scratch = _fill_cells_around_scratch!(
          faces_around_scratch,
          node,
          node_to_faces_data,
          node_to_faces_ptrs)
        _set_intersection!(
          faces_around,faces_around_scratch,
          nfaces_around,nfaces_around_scratch)
      end
    end

    for face in faces_around
      if face != UNSET
        gface_to_face[gface] = face
        break
      end
    end

  end

end
