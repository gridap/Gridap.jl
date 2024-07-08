
const version = "v0.15"

"""
    struct UnstructuredDiscreteModel{Dc,Dp,Tp,B} <: DiscreteModel{Dc,Dp}
      grid::UnstructuredGrid{Dc,Dp,Tp,B}
      grid_topology::UnstructuredGridTopology{Dc,Dp,Tp,B}
      face_labeling::FaceLabeling
    end
"""
struct UnstructuredDiscreteModel{Dc,Dp,Tp,B} <: DiscreteModel{Dc,Dp}
  grid::UnstructuredGrid{Dc,Dp,Tp,B}
  grid_topology::UnstructuredGridTopology{Dc,Dp,Tp,B}
  face_labeling::FaceLabeling
end

"""
    UnstructuredDiscreteModel(grid::Grid)
"""
function UnstructuredDiscreteModel(grid::Grid)
  _grid = UnstructuredGrid(grid)
  topo = UnstructuredGridTopology(_grid)
  nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
  labels = FaceLabeling(nfaces)
  UnstructuredDiscreteModel(_grid,topo,labels)
end

"""
"""
function UnstructuredDiscreteModel(model::DiscreteModel)
  grid = UnstructuredGrid(get_grid(model))
  topo = UnstructuredGridTopology(get_grid_topology(model))
  labels = get_face_labeling(model)
  UnstructuredDiscreteModel(grid,topo,labels)
end

function UnstructuredDiscreteModel(model::UnstructuredDiscreteModel)
  model
end

# Implementation of the interface

get_grid(model::UnstructuredDiscreteModel) = model.grid

get_grid_topology(model::UnstructuredDiscreteModel) = model.grid_topology

get_face_labeling(model::UnstructuredDiscreteModel) = model.face_labeling

# Io

function to_dict(model::UnstructuredDiscreteModel)
  dict = Dict{Symbol,Any}()
  grid = get_grid(model)
  labels = get_face_labeling(model)
  dict[:grid] = to_dict(grid)
  dict[:labeling] = to_dict(labels)
  dict[:vertex_node] = get_vertex_node(model)
  dict[:version] = version
  D = num_cell_dims(model)
  dict[:D] = D
  topo = get_grid_topology(model)
  for d in 1:D
    k = Symbol("face_vertices_$d")
    dict[k] = to_dict(get_faces(topo,d,0))
  end
  dict
end

function check_dict(::Type{UnstructuredDiscreteModel},dict::Dict{Symbol,Any})

  if ! haskey(dict,:version)
    @unreachable "Cannot convert Dict to UnstructuredDiscreteModel for legacy format"
  end

  v = VersionNumber(dict[:version])
  if v < VersionNumber(version)
    @unreachable """\n
    Cannot convert Dict to UnstructuredDiscreteModel for outdated format $(dict[:version]).
    The required format is $version. Regenerate your dictiorary (typically a json file)
    holding the DiscreteModel with the Gridap version you are currently using.
    """
  end

end

function from_dict(T::Type{UnstructuredDiscreteModel},dict::Dict{Symbol,Any})

  check_dict(T,dict)

  grid = from_dict(UnstructuredGrid,dict[:grid])
  labeling = from_dict(FaceLabeling,dict[:labeling])
  D::Int = dict[:D]
  vertex_to_node::Vector{Int32} = dict[:vertex_node]
  orientation = OrientationStyle(grid)
  polytopes = map(get_polytope, get_reffes(grid))
  cell_type = get_cell_type(grid)
  vertex_coordinates = get_node_coordinates(grid)[vertex_to_node]

  nvertices = length(vertex_to_node)
  d_dface_to_vertices = Vector{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,D+1)
  d_dface_to_vertices[0+1] = identity_table(Int32,Int32,nvertices)
  for d in 1:D
    k = Symbol("face_vertices_$d")
    dface_to_vertices = from_dict(Table{Int32,Vector{Int32},Vector{Int32}},dict[k])
    d_dface_to_vertices[d+1] = dface_to_vertices
  end

  topo = UnstructuredGridTopology(
    vertex_coordinates,
    d_dface_to_vertices,
    cell_type,
    polytopes,
    orientation)

  UnstructuredDiscreteModel(grid,topo,labeling)
end

function simplexify(model::UnstructuredDiscreteModel;kwargs...)
  grid = get_grid(model)
  grid_topology = get_grid_topology(model)
  cell_vertices = get_cell_vertices(grid_topology)
  cell_nodes = get_cell_node_ids(grid)
  nnodes = num_nodes(grid)
  nvertices = num_vertices(grid_topology)
  labels = get_face_labeling(model)
  tgrid = simplexify(grid;kwargs...)
  cell_ctype = get_cell_type(grid)
  tcell_nodes = get_cell_node_ids(tgrid)
  tcell_tctype = get_cell_type(tgrid)
  function get_lvertex_lnodes(reffe)
    face_nodes = get_face_own_nodes(reffe)
    poly = get_polytope(reffe)
    first.(face_nodes[get_dimrange(poly,0)])
  end
  reffes = get_reffes(grid)
  ctype_lvertex_lnode = map(get_lvertex_lnodes,reffes)
  treffes = get_reffes(tgrid)
  tctype_ltvertex_ltnode = map(get_lvertex_lnodes,treffes)
  tcell_vertices, vertex_node = _prepare_tvertices(
    cell_vertices,
    cell_nodes,
    nvertices,
    nnodes,
    ctype_lvertex_lnode,
    cell_ctype,
    tcell_nodes,
    tctype_ltvertex_ltnode,
    tcell_tctype)
  tgrid_topology = GridTopology(tgrid, tcell_vertices, vertex_node)
  tfacelabels = _generate_tfacelabels(grid_topology, tgrid_topology, reffes, labels; kwargs... )
  UnstructuredDiscreteModel(tgrid,tgrid_topology,tfacelabels)
end

const UNSET_ID = Int32(0)

function _prepare_tvertices(
  cell_vertices,
  cell_nodes,
  nvertices,
  nnodes,
  ctype_lvertex_lnode,
  cell_ctype,
  tcell_nodes,
  tctype_ltvertex_ltnode,
  tcell_tctype)

  vertex_node = fill(UNSET_ID,nvertices)
  node_vertex = fill(UNSET_ID,nnodes)
  c1 = array_cache(cell_vertices)
  c2 = array_cache(cell_nodes)
  ncells = length(cell_nodes)
  for cell in 1:ncells
    ctype = cell_ctype[cell]
    lvertex_lnode = ctype_lvertex_lnode[ctype]
    vertices = getindex!(c1,cell_vertices,cell)
    nodes = getindex!(c2,cell_nodes,cell)
    for (lvertex,lnode) in enumerate(lvertex_lnode)
      v = vertices[lvertex]
      n = nodes[lnode]
      #if vertex_node[v] != UNSET_ID
      vertex_node[v] = n
      #end
      node_vertex[n] = v
    end
  end

  ntcells = length(tcell_nodes)
  tcell_vertices_ptrs = zeros(Int32,ntcells+1)
  for tcell in 1:ntcells
    tctype = tcell_tctype[tcell]
    ltvertex_ltnode = tctype_ltvertex_ltnode[tctype]
    tcell_vertices_ptrs[tcell+1] = length(ltvertex_ltnode)
  end
  length_to_ptrs!(tcell_vertices_ptrs)
  ndata = tcell_vertices_ptrs[end]-1
  tcell_vertices_data = zeros(Int32,ndata)
  c3 = array_cache(tcell_nodes)
  for tcell in 1:ntcells
    nodes = getindex!(c3,tcell_nodes,tcell)
    tctype = tcell_tctype[tcell]
    ltvertex_ltnode = tctype_ltvertex_ltnode[tctype]
    pini = tcell_vertices_ptrs[tcell]-1
    for (ltvertex,ltnode) in enumerate(ltvertex_ltnode)
      n = nodes[ltnode]
      v = node_vertex[n]
      p = pini + ltvertex
      tcell_vertices_data[p] = v
    end
  end

  tcell_vertices = Table(tcell_vertices_data,tcell_vertices_ptrs)
  tcell_vertices, vertex_node
end

function _generate_tfacelabels(grid_topology, tgrid_topology, reffes, facelabels; kwargs... )

  @notimplementedif length(reffes) != 1
  reffe = first(reffes)
  p = get_polytope(reffe)
  ltcell_to_lnodes, simplex = simplexify(p;kwargs...)

  D = num_cell_dims(grid_topology)

  dim_to_tface_to_label = [ fill(UNSET_ID,num_faces(tgrid_topology,d)) for d in 0:D ]
  dim_to_face_to_label = facelabels.d_to_dface_to_entity

  _fill_dim_to_tface_to_label!(
    dim_to_tface_to_label,
    dim_to_face_to_label,
    tgrid_topology,
    grid_topology,
    ltcell_to_lnodes)

  FaceLabeling(
    dim_to_tface_to_label, facelabels.tag_to_entities, facelabels.tag_to_name)

end

function  _fill_dim_to_tface_to_label!(
  dim_to_tface_to_label,
  dim_to_face_to_label,
  tgrid_topology,
  grid_topology,
  ltcell_to_lnodes)

  tpolytope = first(get_polytopes(tgrid_topology))
  polytope = first(get_polytopes(grid_topology))

  D = length(dim_to_face_to_label)-1
  d = 0
  dim_to_tface_to_label[d+1] = dim_to_face_to_label[d+1]

  for d in 1:(D-1)

    cell_to_faces = get_faces(grid_topology,D,d)
    tcell_to_tfaces = get_faces(tgrid_topology,D,d)

    ntfaces = maximum(tcell_to_tfaces.data)

    ltface_to_ltnodes = get_faces(tpolytope,d,0)
    lface_to_lnodes = get_faces(polytope,d,0)

    tface_to_face = _generate_tface_to_face(
      cell_to_faces.data,
      cell_to_faces.ptrs,
      tcell_to_tfaces.data,
      tcell_to_tfaces.ptrs,
      ltcell_to_lnodes,
      ltface_to_ltnodes,
      lface_to_lnodes,
      ntfaces)

    _update_labels!(
      dim_to_tface_to_label[d+1],dim_to_face_to_label[d+1],tface_to_face)

  end

  d = D
  ncells = length(dim_to_face_to_label[d+1])
  nltcells = length(ltcell_to_lnodes)
  tcell_to_cell = _generate_tcell_to_cell(ncells,nltcells)
  _update_labels!(
    dim_to_tface_to_label[d+1],dim_to_face_to_label[d+1],tcell_to_cell)

  for d = 1:(D-1)

    for j in (d+1):D

      dface_to_jfaces = get_faces(tgrid_topology,d,j)
      dface_to_label = dim_to_tface_to_label[d+1]
      jface_to_label = dim_to_tface_to_label[j+1]
      cache = array_cache(dface_to_jfaces)
      _fix_dface_to_label!(dface_to_label,jface_to_label,dface_to_jfaces,cache)

    end

  end

end

function _fix_dface_to_label!(dface_to_label,jface_to_label,dface_to_jfaces,cache)

  ndfaces = length(dface_to_label)
  @assert ndfaces == length(dface_to_jfaces)

  for dface in 1:ndfaces

    dlabel = dface_to_label[dface]
    if dlabel != UNSET_ID
      continue
    end

    jfaces = getindex!(cache,dface_to_jfaces,dface)
    for jface in jfaces
      jlabel = jface_to_label[jface]
      if jlabel != UNSET_ID
        dface_to_label[dface] = jlabel
        break
      end
    end

  end

end

function _update_labels!(tface_to_label,face_to_label,tface_to_face)
  for tface in 1:length(tface_to_label)
    face = tface_to_face[tface]
    if face != 0
      tface_to_label[tface] = face_to_label[face]
    end
  end
end

function _generate_tcell_to_cell(ncells,nltcells)
  ntcells = ncells * nltcells
  tcell_to_cell = Vector{Int32}(undef,ntcells)
  tcell = 1
  for cell in 1:ncells
    for ltcell in 1:nltcells
      tcell_to_cell[tcell] = cell
      tcell += 1
    end
  end
  tcell_to_cell
end

