
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
    holding the DiscreteModel with the Gridap version you are currenlty using.
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

function simplexify(model::UnstructuredDiscreteModel)

  grid = get_grid(model)
  grid_topology = get_grid_topology(model)
  labels = get_face_labeling(model)

  tgrid = simplexify(grid)
  tgrid_topology = GridTopology(tgrid)

  reffes = get_reffes(grid)
  tfacelabels = _generate_tfacelabels(grid_topology, tgrid_topology, reffes, labels)

  UnstructuredDiscreteModel(tgrid,tgrid_topology,tfacelabels)

end

const UNSET_ID = Int32(0)

function _generate_tfacelabels(grid_topology, tgrid_topology, reffes, facelabels)

  @notimplementedif length(reffes) != 1
  reffe = first(reffes)
  p = get_polytope(reffe)
  ltcell_to_lnodes, simplex = simplexify(p)

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

