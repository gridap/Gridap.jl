
"""
    struct UnstructuredDiscreteModel{Dc,Dp} <: DiscreteModel{Dc,Dp}
      # private fields
    end

Discrete model for any type of unstructured discretization
"""
struct UnstructuredDiscreteModel{Dc,Dp} <: DiscreteModel{Dc,Dp}
  num_nodes::Int
  vertex_to_node::Vector{Int}
  node_to_face_owner::Vector{Int}
  d_to_num_dfaces::Vector{Int}
  d_to_dface_to_nodes::Vector{Table{Int,Int32}}
  d_to_dface_to_isboundary::Vector{Vector{Bool}}
  d_to_dface_to_reffe_type::Vector{Vector{Int8}}
  d_to_dface_to_polytope_type::Vector{Vector{Int8}}
  n_m_to_nface_to_mfaces::Matrix{Table{Int,Int32}}
  d_to_reffes::Vector{Vector{NodalReferenceFE}}
  d_to_polytopes::Vector{Vector{Polytope}}
  node_coordinates::Vector{Point{Dp,Float64}}

  function UnstructuredDiscreteModel(grid::UnstructuredGrid)
    Dc = num_cell_dims(grid)
    Dp = num_point_dims(grid)
    fields = _init_fields(grid)

    ( nnodes,
      vertex_to_node,
      node_to_face_owner,
      d_to_num_dfaces,
      d_to_dface_to_nodes,
      d_to_dface_to_isboundary,
      d_to_dface_to_reffe_type,
      d_to_dface_to_polytope_type,
      n_m_to_nface_to_mfaces,
      d_to_reffes,
      d_to_polytopes,
      node_coordinates) = fields

    new{Dc,Dp}(
     nnodes,
     vertex_to_node,
     node_to_face_owner,
     d_to_num_dfaces,
     d_to_dface_to_nodes,
     d_to_dface_to_isboundary,
     d_to_dface_to_reffe_type,
     d_to_dface_to_polytope_type,
     n_m_to_nface_to_mfaces,
     d_to_reffes,
     d_to_polytopes,
     node_coordinates)
  end
end

"""
    UnstructuredDiscreteModel(trian::ConformingTriangulation)
"""
function UnstructuredDiscreteModel(trian::ConformingTriangulation)
  grid = UnstructuredGrid(trian)
  UnstructuredDiscreteModel(grid)
end

function _init_fields(grid)

  out = generate_cell_to_vertices(grid)
  cell_to_vertices, vertex_to_node, node_to_face_owner = out
  cell_to_type = get_cell_type(grid)
  d = num_cell_dims(grid)
  n = d+1
  nnodes = length(node_to_face_owner)
  nvertices = length(vertex_to_node)
  ncells = num_cells(grid)
  vertex_to_cells = generate_cells_around(cell_to_vertices,nvertices)
  node_coordinates = get_node_coordinates(grid)

  d_to_num_dfaces = fill(UNSET,n)
  d_to_dface_to_nodes = Vector{Table{Int,Int32}}(undef,n)
  d_to_dface_to_isboundary = Vector{Vector{Bool}}(undef,n)
  d_to_dface_to_reffe_type = Vector{Vector{Int8}}(undef,n)
  d_to_dface_to_polytope_type = Vector{Vector{Int8}}(undef,n)
  n_m_to_nface_to_mfaces = Matrix{Table{Int,Int32}}(undef,n,n)
  d_to_reffes = Vector{Vector{NodalReferenceFE}}(undef,n)
  d_to_polytopes = Vector{Vector{Polytope}}(undef,n)

  n_m_to_nface_to_mfaces[0+1,0+1] = identity_table(Int,Int32,nvertices)
  n_m_to_nface_to_mfaces[0+1,d+1] = vertex_to_cells
  n_m_to_nface_to_mfaces[d+1,0+1] = cell_to_vertices
  n_m_to_nface_to_mfaces[d+1,d+1] = identity_table(Int,Int32,ncells)
  d_to_num_dfaces[0+1] = length(vertex_to_node)
  d_to_num_dfaces[d+1] = length(cell_to_vertices)
  d_to_dface_to_nodes[d+1] = get_cell_nodes(grid)
  d_to_dface_to_reffe_type[d+1] = cell_to_type
  d_to_dface_to_polytope_type[d+1] = cell_to_type
  d_to_reffes[d+1] = get_reffes(grid)
  d_to_polytopes[d+1] = map(get_polytope,get_reffes(grid))

  fields = (
    nnodes,
    vertex_to_node,
    node_to_face_owner,
    d_to_num_dfaces,
    d_to_dface_to_nodes,
    d_to_dface_to_isboundary,
    d_to_dface_to_reffe_type,
    d_to_dface_to_polytope_type,
    n_m_to_nface_to_mfaces,
    d_to_reffes,
    d_to_polytopes,
    node_coordinates)

  fields
end


function num_faces(g::UnstructuredDiscreteModel,d::Integer)
  _generate_cell_to_faces!(g,d)
  g.d_to_num_dfaces[d+1]
end

function num_nodes(g::UnstructuredDiscreteModel)
  g.num_nodes
end

function get_faces(g::UnstructuredDiscreteModel,dimfrom::Integer,dimto::Integer)
  D = num_cell_dims(g)

  if dimfrom==0
    if dimto==0 || dimto==D
      nothing
    else
      _generate_face_to_vertices!(g,dimto)
    end

  elseif dimfrom == D
    if dimto==0 || dimto==D
      nothing
    else
      _generate_cell_to_faces!(g,dimto)
    end

  else

    if dimto==0
      _generate_face_to_vertices!(g,dimfrom)
    elseif dimto==D
      _generate_cell_to_faces!(g,dimfrom)
    elseif dimto==dimfrom
      _generate_face_to_face!(g,dimfrom)
    elseif dimfrom > dimto
      _generate_nface_to_mface!(g,dimfrom,dimto)
    else
      _generate_nface_to_mface!(g,dimto,dimfrom)
    end

  end

  g.n_m_to_nface_to_mfaces[dimfrom+1,dimto+1]
end

function _generate_cell_to_faces!(model,dimto)

  D = num_cell_dims(model)

  if isassigned(model.n_m_to_nface_to_mfaces,D+1,dimto+1)
    return
  end
  if isassigned(model.n_m_to_nface_to_mfaces,dimto+1,0+1)
    @notimplemented

  else

    cell_to_vertices = model.n_m_to_nface_to_mfaces[D+1,0+1]
    vertex_to_cells = model.n_m_to_nface_to_mfaces[0+1,D+1]
    cell_to_cell_type = model.d_to_dface_to_polytope_type[D+1]
    polytopes = model.d_to_polytopes[D+1]
    cell_type_to_lface_to_lvertices = map( (p)->get_faces(p,dimto,0), polytopes )

    cell_to_faces = generate_cell_to_faces(
        cell_to_vertices,
        cell_type_to_lface_to_lvertices,
        cell_to_cell_type,
        vertex_to_cells)

    faces_to_cells = generate_cells_around(cell_to_faces)
    nfaces = length(faces_to_cells)

    model.n_m_to_nface_to_mfaces[D+1,dimto+1] = cell_to_faces
    model.n_m_to_nface_to_mfaces[dimto+1,D+1] = faces_to_cells
    model.d_to_num_dfaces[dimto+1] = nfaces

  end

  nothing

end

function _generate_face_to_vertices!(model,dimfrom)

  if isassigned(model.n_m_to_nface_to_mfaces,dimfrom+1,0+1)
    return
  end

  D = num_cell_dims(model)

  cell_to_vertices = get_faces(model,D,0)
  cell_to_faces = get_faces(model,D,dimfrom)
  cell_to_ctype = model.d_to_dface_to_polytope_type[D+1]
  polytopes = model.d_to_polytopes[D+1]
  ctype_to_lface_to_lvertices = map( (p)->get_faces(p,dimfrom,0), polytopes )

  nfaces = model.d_to_num_dfaces[dimfrom+1]

  face_to_vertices = generate_face_to_vertices(
    cell_to_vertices,
    cell_to_faces,
    cell_to_ctype,
    ctype_to_lface_to_lvertices,
    nfaces)

  nvertices = length(model.vertex_to_node)

  vertex_to_faces = generate_cells_around(face_to_vertices,nvertices)

  model.n_m_to_nface_to_mfaces[dimfrom+1,0+1] = face_to_vertices
  model.n_m_to_nface_to_mfaces[0+1,dimfrom+1] = vertex_to_faces

  nothing

end

function _generate_face_to_face!(model,d)

  if isassigned(model.n_m_to_nface_to_mfaces,d+1,d+1)
    return
  end

  nfaces = num_faces(model,d)
  id = identity_table(Int,Int32,nfaces)
  model.n_m_to_nface_to_mfaces[d+1,d+1] = id

  return

end

function _generate_nface_to_mface!(model,n,m)
  @assert n > m

  if isassigned(model.n_m_to_nface_to_mfaces,n+1,m+1)
    return
  end

  nface_to_vertices = get_faces(model,n,0)
  vertex_to_mfaces = get_faces(model,0,m)
  nface_to_nftype = get_face_polytope_type(model,n)

  polytopes = get_polytopes(Polytope{n},model)
  nftype_to_lmface_to_lvertices = map( (p)->get_faces(p,m,0), polytopes)

  nface_to_mfaces = find_cell_to_faces(
    nface_to_vertices,
    nftype_to_lmface_to_lvertices,
    nface_to_nftype,
    vertex_to_mfaces)

  num_mfaces = num_faces(model,m)
  mface_to_nfaces = generate_cells_around(nface_to_mfaces,num_mfaces)

  model.n_m_to_nface_to_mfaces[n+1,m+1] = nface_to_mfaces
  model.n_m_to_nface_to_mfaces[m+1,n+1] = mface_to_nfaces

  nothing

end

function get_vertex_node(g::UnstructuredDiscreteModel)
  g.vertex_to_node
end

function get_node_face_owner(g::UnstructuredDiscreteModel)
  g.node_to_face_owner
end

function get_face_nodes(g::UnstructuredDiscreteModel,d::Integer)
  _generate_face_to_nodes!(g,d)
  g.d_to_dface_to_nodes[d+1]
end

function get_cell_nodes(g::UnstructuredDiscreteModel)
  D = num_cell_dims(g)
  g.d_to_dface_to_nodes[D+1]
end

function _generate_face_to_nodes!(model,dimfrom)

  if isassigned(model.d_to_dface_to_nodes,dimfrom+1)
    return
  end

  D = num_cell_dims(model)

  cell_to_nodes = get_cell_nodes(model)
  cell_to_faces = get_faces(model,D,dimfrom)
  cell_to_ctype = model.d_to_dface_to_reffe_type[D+1]
  reffes = model.d_to_reffes[D+1]
  ctype_to_lface_to_lnodes = map( (p)->get_face_nodes(p)[get_dimranges(get_polytope(p))[dimfrom+1]], reffes )

  nfaces = num_faces(model,dimfrom)

  face_to_nodes = generate_face_to_vertices(
    cell_to_nodes,
    cell_to_faces,
    cell_to_ctype,
    ctype_to_lface_to_lnodes,
    nfaces)

  model.d_to_dface_to_nodes[dimfrom+1] = face_to_nodes

  nothing

end

function get_isboundary_face(g::UnstructuredDiscreteModel,d::Integer)
  D = num_cell_dims(g)
  if d == D-1
    _generate_isboundary_facet!(g)
  else
    _generate_isboundary_face!(g,d)
  end
  g.d_to_dface_to_isboundary[d+1]
end

function _generate_isboundary_facet!(model)

  D = num_cell_dims(model)
  d = D-1
  if isassigned(model.d_to_dface_to_isboundary,d+1)
    return
  end

  face_to_cells = get_faces(model,d,D)
  face_to_isboundary = generate_facet_to_isboundary(face_to_cells)
  model.d_to_dface_to_isboundary[d+1] = face_to_isboundary

  nothing

end

function _generate_isboundary_face!(model,d)

  if isassigned(model.d_to_dface_to_isboundary,d+1)
    return
  end

  D = num_cell_dims(model)

  if d==D
    _generate_cell_to_faces!(model,D-1)
  elseif d==0
    _generate_face_to_vertices!(model,D-1)
  end

  if !isassigned(model.n_m_to_nface_to_mfaces,d+1,D-1+1)

    cell_to_faces = get_faces(model,D,d)
    cell_to_facets = get_faces(model,D,D-1)
    cell_to_ctype = get_face_polytope_type(model,D)
    polytopes = get_polytopes(Polytope{D},model)
    ctype_to_lface_to_lfacets = map( (p)->get_faces(p,d,D-1), polytopes)
    facet_to_isboundary = get_isboundary_face(model,D-1)
    nfaces = num_faces(model,d)

    face_to_isboundary = generate_face_to_isboundary_from_cells(
      facet_to_isboundary,
      cell_to_faces,
      cell_to_facets,
      cell_to_ctype,
      ctype_to_lface_to_lfacets,
      nfaces)

  else
    facet_to_isboundary = get_isboundary_face(model,D-1)
    face_to_facets = get_faces(model,d,D-1)
    face_to_isboundary = generate_face_to_isboundary(
      facet_to_isboundary, face_to_facets)
  end

  model.d_to_dface_to_isboundary[d+1] = face_to_isboundary

  nothing

end

function get_face_reffe_type(g::UnstructuredDiscreteModel,d::Integer)
  _generate_face_reffes!(g,d)
  g.d_to_dface_to_reffe_type[d+1]
end

function get_face_polytope_type(g::UnstructuredDiscreteModel,d::Integer)
  _generate_face_polytopes!(g,d)
  g.d_to_dface_to_polytope_type[d+1]
end

function get_reffes(T::Type{<:ReferenceFE{d}},g::UnstructuredDiscreteModel) where d
  _generate_face_reffes!(g,d)
  reffes::Vector{T} = g.d_to_reffes[d+1]
  reffes
end

function get_polytopes(T::Type{<:Polytope{d}},g::UnstructuredDiscreteModel) where d
  _generate_face_polytopes!(g,d)
  polytopes::Vector{T} = g.d_to_polytopes[d+1]
  polytopes
end

function get_node_coordinates(g::UnstructuredDiscreteModel)
  g.node_coordinates
end

function _generate_face_reffes!(model,d)

  if isassigned(model.d_to_dface_to_reffe_type,d+1)
    return
  end

  D = num_cell_dims(model)

  ctype_to_reffe = model.d_to_reffes[D+1]

  ctype_to_lftype_to_refface = [ get_reffes(ReferenceFE{d},reffe) for reffe in ctype_to_reffe]
  ctype_to_lface_to_lftype = [ get_face_type(reffe,d) for reffe in ctype_to_reffe]

  t = _generate_ftype_to_refface(Val{d}(),ctype_to_lftype_to_refface,ctype_to_lface_to_lftype)
  ftype_to_refface, ctype_to_lface_to_ftype = t

  cell_to_faces = get_faces(model,D,d)
  cell_to_ctype = model.d_to_dface_to_reffe_type[D+1]

  nfaces = num_faces(model,d)

  face_to_ftype = generate_face_to_face_type(
    cell_to_faces, cell_to_ctype, ctype_to_lface_to_ftype, nfaces)

  model.d_to_dface_to_reffe_type[d+1] = face_to_ftype
  model.d_to_reffes[d+1] = ftype_to_refface

  nothing

end

function _generate_face_polytopes!(model,d)

  if isassigned(model.d_to_dface_to_polytope_type,d+1)
    return
  end

  D = num_cell_dims(model)

  ctype_to_polytope = model.d_to_polytopes[D+1]

  ctype_to_lftype_to_refface = [ get_reffaces(Polytope{d},polytope) for polytope in ctype_to_polytope]
  ctype_to_lface_to_lftype = [ get_face_type(polytope,d) for polytope in ctype_to_polytope]

  t = _generate_ftype_to_refface(Val{d}(),ctype_to_lftype_to_refface,ctype_to_lface_to_lftype)
  ftype_to_refface, ctype_to_lface_to_ftype = t

  cell_to_faces = get_faces(model,D,d)
  cell_to_ctype = model.d_to_dface_to_polytope_type[D+1]

  nfaces = num_faces(model,d)

  face_to_ftype = generate_face_to_face_type(
    cell_to_faces, cell_to_ctype, ctype_to_lface_to_ftype, nfaces)

  model.d_to_dface_to_polytope_type[d+1] = face_to_ftype
  model.d_to_polytopes[d+1] = ftype_to_refface

  nothing

end

# Low level helpers

include("GridOperations.jl")

