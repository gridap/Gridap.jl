
"""
    struct UnstructuredDiscreteModel{Dc,Dp,B} <: DiscreteModel{Dc,Dp}
      # private fields
    end

Discrete model for any type of unstructured discretization
"""
struct UnstructuredDiscreteModel{Dc,Dp,B} <: DiscreteModel{Dc,Dp}
  num_nodes::Int
  vertex_to_node::Vector{Int}
  node_to_face_owner::Vector{Int}
  d_to_num_dfaces::Vector{Int}
  d_to_dface_to_nodes::Vector{Table{Int,Int32}}
  d_to_dface_to_own_nodes::Vector{Table{Int,Int32}}
  d_to_dface_to_isboundary::Vector{Vector{Bool}}
  d_to_dface_to_reffe_type::Vector{Vector{Int8}}
  d_to_dface_to_polytope_type::Vector{Vector{Int8}}
  n_m_to_nface_to_mfaces::Matrix{Table{Int,Int32}}
  d_to_reffes::Vector{Vector{NodalReferenceFE}}
  d_to_polytopes::Vector{Vector{Polytope}}
  node_coordinates::Vector{Point{Dp,Float64}}
  labels::FaceLabeling
end

function UnstructuredDiscreteModel(grid::UnstructuredGrid)
  Dc = num_cell_dims(grid)
  Dp = num_point_dims(grid)
  fields = _init_fields(grid)
  B = is_oriented(grid)

  ( nnodes,
    vertex_to_node,
    node_to_face_owner,
    d_to_num_dfaces,
    d_to_dface_to_nodes,
    d_to_dface_to_own_nodes,
    d_to_dface_to_isboundary,
    d_to_dface_to_reffe_type,
    d_to_dface_to_polytope_type,
    n_m_to_nface_to_mfaces,
    d_to_reffes,
    d_to_polytopes,
    node_coordinates,
    labels) = fields

  UnstructuredDiscreteModel{Dc,Dp,B}(
   nnodes,
   vertex_to_node,
   node_to_face_owner,
   d_to_num_dfaces,
   d_to_dface_to_nodes,
   d_to_dface_to_own_nodes,
   d_to_dface_to_isboundary,
   d_to_dface_to_reffe_type,
   d_to_dface_to_polytope_type,
   n_m_to_nface_to_mfaces,
   d_to_reffes,
   d_to_polytopes,
   node_coordinates,
   labels)
end

"""
    UnstructuredDiscreteModel(trian::ConformingTriangulation)
"""
function UnstructuredDiscreteModel(trian::ConformingTriangulation)
  grid = UnstructuredGrid(trian)
  UnstructuredDiscreteModel(grid)
end

OrientationStyle(
  ::Type{UnstructuredDiscreteModel{Dc,Dp,B}}) where {Dc,Dp,B} = Val{B}()

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
  d_to_dface_to_own_nodes = Vector{Table{Int,Int32}}(undef,n)
  d_to_dface_to_isboundary = Vector{Vector{Bool}}(undef,n)
  d_to_dface_to_reffe_type = Vector{Vector{Int8}}(undef,n)
  d_to_dface_to_polytope_type = Vector{Vector{Int8}}(undef,n)
  n_m_to_nface_to_mfaces = Matrix{Table{Int,Int32}}(undef,n,n)
  d_to_reffes = Vector{Vector{NodalReferenceFE}}(undef,n)
  d_to_polytopes = Vector{Vector{Polytope}}(undef,n)
  d_to_dface_to_entity = Vector{Vector{Int32}}(undef,n)
  tag_to_entities = Vector{Int32}[]
  tag_to_name = String[]
  labels = FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)

  n_m_to_nface_to_mfaces[0+1,0+1] = identity_table(Int,Int32,nvertices)
  n_m_to_nface_to_mfaces[0+1,d+1] = vertex_to_cells
  n_m_to_nface_to_mfaces[d+1,0+1] = cell_to_vertices
  n_m_to_nface_to_mfaces[d+1,d+1] = identity_table(Int,Int32,ncells)
  d_to_num_dfaces[0+1] = length(vertex_to_node)
  d_to_num_dfaces[d+1] = length(cell_to_vertices)
  d_to_dface_to_nodes[d+1] = get_cell_nodes(grid)
  d_to_dface_to_own_nodes[0+1] = Table(vertex_to_node,collect(Int32,1:(nvertices+1)))
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
    d_to_dface_to_own_nodes,
    d_to_dface_to_isboundary,
    d_to_dface_to_reffe_type,
    d_to_dface_to_polytope_type,
    n_m_to_nface_to_mfaces,
    d_to_reffes,
    d_to_polytopes,
    node_coordinates,
    labels)

  fields
end

"""
"""
function UnstructuredDiscreteModel(g::DiscreteModel)

  Dc = num_cell_dims(g)
  Dp = num_point_dims(g)
  B = is_oriented(g)
  D = Dc

  nnodes = num_nodes(g)
  vertex_to_node = get_vertex_node(g)
  node_to_face_owner = get_node_face_owner(g)
  d_to_num_dfaces = [ num_faces(g,d) for d in 0:D ]
  d_to_dface_to_nodes = [ Table(get_face_nodes(g,d)) for d in 0:D ]
  d_to_dface_to_own_nodes = [ Table(get_face_own_nodes(g,d)) for d in 0:D ]
  d_to_dface_to_isboundary = [ get_isboundary_face(g,d) for d in 0:D ]
  d_to_dface_to_reffe_type = [ get_face_reffe_type(g,d) for d in 0:D ]
  d_to_dface_to_polytope_type = [ get_face_polytope_type(g,d) for d in 0:D ]
  n_m_to_nface_to_mfaces = [ get_faces(g,i,j) for i in 0:D, j in 0:D ]
  d_to_reffes = [get_reffes(NodalReferenceFE{d},g) for d in 0:D]
  d_to_polytopes = [get_polytopes(Polytope{d},g) for d in 0:D]
  node_coordinates = get_node_coordinates(g)
  labels = get_face_labeling(g)

  UnstructuredDiscreteModel{Dc,Dp,B}(
   nnodes,
   vertex_to_node,
   node_to_face_owner,
   d_to_num_dfaces,
   d_to_dface_to_nodes,
   d_to_dface_to_own_nodes,
   d_to_dface_to_isboundary,
   d_to_dface_to_reffe_type,
   d_to_dface_to_polytope_type,
   n_m_to_nface_to_mfaces,
   d_to_reffes,
   d_to_polytopes,
   node_coordinates,
   labels)

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

function get_face_own_nodes(g::UnstructuredDiscreteModel,d::Integer)
  _generate_face_to_own_nodes!(g,d)
  g.d_to_dface_to_own_nodes[d+1]
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

function _generate_face_to_own_nodes!(model,dimfrom)

  if isassigned(model.d_to_dface_to_own_nodes,dimfrom+1)
    return
  end

  D = num_cell_dims(model)

  cell_to_nodes = get_cell_nodes(model)
  cell_to_faces = get_faces(model,D,dimfrom)
  cell_to_ctype = model.d_to_dface_to_reffe_type[D+1]
  reffes = model.d_to_reffes[D+1]
  ctype_to_lface_to_lnodes = map( (p)->get_face_own_nodes(p)[get_dimranges(get_polytope(p))[dimfrom+1]], reffes )

  nfaces = num_faces(model,dimfrom)

  face_to_nodes = generate_face_to_vertices(
    cell_to_nodes,
    cell_to_faces,
    cell_to_ctype,
    ctype_to_lface_to_lnodes,
    nfaces)

  model.d_to_dface_to_own_nodes[dimfrom+1] = face_to_nodes

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

function get_face_labeling(g::UnstructuredDiscreteModel)
  _generate_face_labeling!(g)
  g.labels
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

function _generate_face_labeling!(model)

  if isassigned(model.labels.d_to_dface_to_entity,0+1)
    return
  end

  for d in 0:num_cell_dims(model)
    model.labels.d_to_dface_to_entity[d+1] = fill(Int32(UNSET),num_faces(model,d))
  end

  nothing
end

function replace_reffes(
  model::UnstructuredDiscreteModel{Dc,Dp,B},
  reffes::Vector{<:NodalReferenceFE}) where {Dc,Dp,B}

  # TODO It assumes that p-conformity is mantained. Check it?

  face_to_own_nodes, n_nodes = _generate_face_to_own_nodes(
    model,reffes,ConformityStyle(model))

  node_to_face_owner  = flatten_partition(face_to_own_nodes,n_nodes)

  cell_to_nodes = _generate_cell_to_nodes(
    face_to_own_nodes, model, reffes, OrientationStyle(model), ConformityStyle(model))

  node_to_coords = _generate_node_to_coords(
    n_nodes, cell_to_nodes, model, reffes, ConformityStyle(model))

  d = num_cell_dims(model)
  n = d+1
  d_to_num_dfaces = model.d_to_num_dfaces
  d_to_dface_to_nodes = Vector{Table{Int,Int32}}(undef,n)
  d_to_dface_to_nodes[d+1] = cell_to_nodes
  d_to_dface_to_own_nodes = [ face_to_own_nodes[range] for range in get_dimranges(model)]
  vertex_to_node = d_to_dface_to_own_nodes[0+1].data
  d_to_dface_to_isboundary = model.d_to_dface_to_isboundary
  d_to_dface_to_reffe_type = Vector{Vector{Int8}}(undef,n)
  d_to_dface_to_reffe_type[d+1] = get_cell_type(model)
  d_to_dface_to_polytope_type = model.d_to_dface_to_polytope_type
  n_m_to_nface_to_mfaces =  model.n_m_to_nface_to_mfaces
  d_to_reffes = Vector{Vector{NodalReferenceFE}}(undef,n)
  d_to_reffes[d+1] = reffes
  d_to_polytopes = model.d_to_polytopes
  node_coordinates = node_to_coords
  labels = model.labels

  UnstructuredDiscreteModel{Dc,Dp,B}(
    n_nodes,
    vertex_to_node,
    node_to_face_owner,
    d_to_num_dfaces,
    d_to_dface_to_nodes,
    d_to_dface_to_own_nodes,
    d_to_dface_to_isboundary,
    d_to_dface_to_reffe_type,
    d_to_dface_to_polytope_type,
    n_m_to_nface_to_mfaces,
    d_to_reffes,
    d_to_polytopes,
    node_coordinates,
    labels)

end

function _generate_face_to_own_nodes(model,reffes,conf)
  @notimplemented
end

function _generate_face_to_own_nodes(model,reffes,conf::RegularConformity)

  D = num_cell_dims(model)
  cell_to_ctype = get_cell_type(model)

  d_to_ctype_to_ldface_to_own_lnodes = [
    [ get_face_own_nodes(reffe,d) for reffe in reffes] for d in 0:D]

  d_to_dface_to_cells = [ get_faces(model,d,D) for d in 0:D]
  d_to_cell_to_dfaces = [ get_faces(model,D,d) for d in 0:D]

  d_to_offset = get_offsets(model)

  n_faces = num_faces(model)

  face_to_own_nodes, n_nodes =  _generate_face_to_own_dofs(
    n_faces,
    cell_to_ctype,
    d_to_cell_to_dfaces,
    d_to_dface_to_cells,
    d_to_offset,
    d_to_ctype_to_ldface_to_own_lnodes)

  (face_to_own_nodes, n_nodes)

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

  for d in 0:D
    cell_to_dfaces = d_to_cell_to_dfaces[d+1]
    dface_to_cells = d_to_dface_to_cells[d+1]
    offset = d_to_offset[d+1]
    ctype_to_ldface_to_own_ldofs = d_to_ctype_to_ldface_to_own_ldofs[d+1]
    ctype_to_ldface_to_num_own_ldofs = map( (x) -> length.(x) ,ctype_to_ldface_to_own_ldofs)
    icell = 1
    dface_to_cell_owner = get_local_item(dface_to_cells,icell)
    dface_to_ldface = find_local_index(dface_to_cell_owner,cell_to_dfaces)

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
  (face_to_own_dofs, n_dofs)
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

function _generate_cell_to_nodes(
  face_to_own_nodes, model, reffes, orientation, conformity)
  @notimplemented
end

function _generate_cell_to_nodes(
  face_to_own_nodes, model, reffes, orientation::Val{true}, ::RegularConformity)

  cell_to_faces = get_cell_faces(model)
  cell_to_ctype = get_cell_type(model)
  ctype_lface_to_own_lnodes = [get_face_own_nodes(reffe) for reffe in reffes]
  ctype_to_num_dofs = map(num_dofs,reffes)

  cell_to_nodes = CellDofsOriented(
    cell_to_faces,
    cell_to_ctype,
    ctype_lface_to_own_lnodes,
    ctype_to_num_dofs,
    face_to_own_nodes)

  Table(cell_to_nodes)

end

function _generate_cell_to_nodes(
  face_to_own_nodes, model, reffes, orientation::Val{false}, ::RegularConformity)

  cell_to_faces = get_cell_faces(model)
  cell_to_lface_to_pindex = get_cell_perm_indices(model)
  cell_to_ctype = get_cell_type(model)
  ctype_to_lface_to_own_lnodes = map( get_face_own_nodes , reffes)
  ctype_to_num_nodes = map( num_nodes , reffes)
  facereffes, face_to_ftype = extract_face_reffes(model,reffes)
  ftype_to_pindex_to_pnodes = map(get_own_dofs_permutations,facereffes)

  cell_to_nodes = CellDofsNonOriented(
    cell_to_faces,
    cell_to_lface_to_pindex,
    cell_to_ctype,
    ctype_to_lface_to_own_lnodes,
    ctype_to_num_nodes,
    face_to_own_nodes,
    face_to_ftype,
    ftype_to_pindex_to_pnodes)

  Table(cell_to_nodes)

end

struct CellDofsOriented{T,P,V<:AbstractVector} <: AbstractVector{Vector{T}}
  cell_to_faces::Table{T,P}
  cell_to_ctype::V
  ctype_to_lface_to_own_ldofs::Vector{Vector{Vector{Int}}}
  ctype_to_num_dofs::Vector{Int}
  face_to_own_dofs::Table{T,P}
end

Base.size(a::CellDofsOriented) = size(a.cell_to_faces)

Base.IndexStyle(::Type{<:CellDofsOriented}) = IndexStyle(Table)

function array_cache(a::CellDofsOriented)
  n_dofs = testitem(a.ctype_to_num_dofs)
  T = eltype(eltype(a))
  v = zeros(T,n_dofs)
  CachedArray(v)
end

function getindex!(cache,a::CellDofsOriented,cell::Integer)
  ctype = a.cell_to_ctype[cell]
  n_dofs = a.ctype_to_num_dofs[ctype]
  setsize!(cache,(n_dofs,))
  dofs = cache.array
  lface_to_own_ldofs = a.ctype_to_lface_to_own_ldofs[ctype]
  p = a.cell_to_faces.ptrs[cell]-1
  for (lface, own_ldofs) in enumerate(lface_to_own_ldofs)
    face = a.cell_to_faces.data[p+lface]
    q = a.face_to_own_dofs.ptrs[face]-1
    for (i,ldof) in enumerate(own_ldofs)
      dof = a.face_to_own_dofs.data[q+i]
      dofs[ldof] = dof
    end
  end
  dofs
end

function Base.getindex(a::CellDofsOriented,cell::Integer)
  cache = array_cache(a)
  getindex!(cache,a,cell)
end

struct CellDofsNonOriented{T,P,C<:AbstractVector,F<:AbstractVector} <:AbstractVector{Vector{T}}
  cell_to_faces::Table{T,P}
  cell_to_lface_to_pindex::Table{T,P}
  cell_to_ctype::C
  ctype_to_lface_to_own_ldofs::Vector{Vector{Vector{Int}}}
  ctype_to_num_dofs::Vector{Int}
  face_to_own_dofs::Table{T,P}
  face_to_ftype::F
  ftype_to_pindex_to_pdofs::Vector{Vector{Vector{Int}}}
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
    ftype = a.face_to_ftype[face]
    pindex = a.cell_to_lface_to_pindex.data[p+lface]
    pdofs = a.ftype_to_pindex_to_pdofs[ftype][pindex]

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

function _generate_node_to_coords(n_nodes, cell_to_lnode_to_node, model,reffes,conformity)
  @notimplemented
end

function _generate_node_to_coords(n_nodes, cell_to_lnode_to_node, model, reffes, ::RegularConformity)

  cell_to_map = get_cell_map(model)
  ctype_to_lnode_to_lcoords = map(get_node_coordinates,reffes)
  cell_to_ctype = get_cell_type(model)
  cell_to_lnode_to_lcoords = _get_cell_data(ctype_to_lnode_to_lcoords, cell_to_ctype)
  cell_to_lnode_to_coords = evaluate(cell_to_map,cell_to_lnode_to_lcoords)

  T = eltype(eltype(get_node_coordinates(model)))
  D = num_point_dims(model)
  node_to_coords = zeros(Point{D,T},n_nodes)

  cache = array_cache(cell_to_lnode_to_coords)

  _fill_node_to_coords!(
    node_to_coords,
    cell_to_lnode_to_node,
    cell_to_lnode_to_coords,
    cache)

  node_to_coords

end

function  _fill_node_to_coords!(
  node_to_coords,
  cell_to_lnode_to_node::Table,
  cell_to_lnode_to_coords,
  cache)
  for cell in 1:length(cell_to_lnode_to_node)
    p = cell_to_lnode_to_node.ptrs[cell]-1
    lnode_to_coords = getindex!(cache,cell_to_lnode_to_coords,cell)
    for (lnode, coords) in enumerate(lnode_to_coords)
      node = cell_to_lnode_to_node.data[p+lnode]
      node_to_coords[node] = coords
    end
  end
end

# Low level helpers

#include("GridOperations.jl")

