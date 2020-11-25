
"""
    struct CartesianDiscreteModel{D,T,F} <: DiscreteModel{D,D}
      # Private Fields
    end
"""
struct CartesianDiscreteModel{D,T,F} <: DiscreteModel{D,D}
  grid::CartesianGrid{D,T,F}
  grid_topology::UnstructuredGridTopology{D,D,T,Oriented}
  face_labeling::FaceLabeling
  @doc """
      CartesianDiscreteModel(desc::CartesianDescriptor)

  Inner constructor
  """
  function CartesianDiscreteModel(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
    grid = CartesianGrid(desc)
    _grid = UnstructuredGrid(grid)
    if any(desc.isperiodic)
      topo = _cartesian_grid_topology_with_periodic_bcs(_grid, desc.isperiodic, desc.partition)
    else
      topo = UnstructuredGridTopology(_grid)
    end
    nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
    labels = FaceLabeling(nfaces)
    _fill_cartesian_face_labeling!(labels,topo)
    new{D,T,F}(grid,topo,labels)
  end

  @doc """
      CartesianDiscreteModel(desc::CartesianDescriptor{D,T,F},
                             cmin::CartesianIndex,
                             cmax::CartesianIndex)

      Builds a CartesianDiscreteModel object which represents a subgrid of
      a (larger) grid represented by desc. This subgrid is described by its
      D-dimensional minimum (cmin) and maximum (cmax) CartesianIndex
      identifiers.

  Inner constructor
  """
  function CartesianDiscreteModel(desc::CartesianDescriptor{D,T,F},
                                  cmin::CartesianIndex,
                                  cmax::CartesianIndex) where {D,T,F}

     suborigin = Tuple(desc.origin) .+ (Tuple(cmin) .- 1) .* desc.sizes
     subpartition = Tuple(cmax) .- Tuple(cmin) .+ 1
     subsizes = desc.sizes
     subdesc =
       CartesianDescriptor(Point(suborigin), subsizes, subpartition; map=desc.map, isperiodic=desc.isperiodic)

     grid = CartesianGrid(subdesc)
     _grid = UnstructuredGrid(grid)
     topo = UnstructuredGridTopology(_grid)
     nfaces = [num_faces(topo, d) for d = 0:num_cell_dims(topo)]
     labels = FaceLabeling(nfaces)
     _fill_subgrid_cartesian_face_labeling!(labels,topo,subdesc,desc,cmin)
     new{D,T,F}(grid, topo, labels)
  end
end

"""
    CartesianDiscreteModel(args...)

Same args needed to construct a `CartesianDescriptor`
"""
function CartesianDiscreteModel(args...; kwargs...)
  desc = CartesianDescriptor(args...; kwargs...)
  CartesianDiscreteModel(desc)
end


"""
    get_cartesian_descriptor(model::CartesianDiscreteModel)
"""
function get_cartesian_descriptor(model::CartesianDiscreteModel)
  get_cartesian_descriptor(model.grid)
end

# Interface

get_grid(model::CartesianDiscreteModel) = model.grid

get_grid_topology(model::CartesianDiscreteModel) = model.grid_topology

get_face_labeling(model::CartesianDiscreteModel) = model.face_labeling

# These needed to be type stable

function get_face_nodes(model::CartesianDiscreteModel,d::Integer)
  face_nodes::Table{Int32,Vector{Int32},Vector{Int32}} = compute_face_nodes(model,d)
  face_nodes
end

function get_face_type(model::CartesianDiscreteModel,d::Integer)
  _, face_to_ftype::Vector{Int8} = compute_reffaces(ReferenceFE{d},model)
  face_to_ftype
end

function get_reffaces(::Type{ReferenceFE{d}},model::CartesianDiscreteModel) where d
  reffaces::Vector{LagrangianRefFE{d}},_ = compute_reffaces(ReferenceFE{d},model)
  reffaces
end

# Helpers

function _fill_cartesian_face_labeling!(labels,topo)
  _fill_cartesian_entities!(labels,topo)
  _add_cartesian_tags!(labels,topo)
end

function _fill_cartesian_entities!(labels,topo)
  D = num_cell_dims(topo)
  d_to_dface_to_entity = labels.d_to_dface_to_entity
  polytope = first(get_polytopes(topo))
  dim_to_offset = get_offsets(polytope)
  interior_id = num_faces(polytope)
  boundary_id = -1
  for d in 0:(D-1)
    face_to_cells = get_faces(topo,d,D)
    cell_to_faces = get_faces(topo,D,d)
    offset = dim_to_offset[d+1]
    _generate_pre_geolabel!(
      d_to_dface_to_entity[d+1],
      face_to_cells,
      cell_to_faces,
      offset,
      interior_id,boundary_id,d,D)
  end
  _fix_geolabels(D, topo, d_to_dface_to_entity, interior_id, boundary_id)
  fill!(d_to_dface_to_entity[end],interior_id)
end

function _fix_geolabels(D, topo, d_to_dface_to_entity, interior_id, boundary_id)
  for d = 0:(D-2)
    for j = (d+1):(D-1)
      dface_to_jfaces = get_faces(topo, d, j)
      dface_to_geolabel = d_to_dface_to_entity[d+1]
      jface_to_geolabel = d_to_dface_to_entity[j+1]
      _fix_dface_geolabels!(
        dface_to_geolabel,
        jface_to_geolabel,
        dface_to_jfaces.data,
        dface_to_jfaces.ptrs,
        interior_id,
        boundary_id,
      )
    end
  end
end

function _add_cartesian_tags!(labels,topo)
  D = num_cell_dims(topo)
  polytope = first(get_polytopes(topo))
  interior_id = num_faces(polytope)
  boundary_ids = collect(1:(interior_id-1))
  for i in boundary_ids
    name = lpad(i,ceil(Int,log10(interior_id)),'0')
    add_tag!(labels,"tag_$(name)",[i])
  end
  add_tag!(labels,"interior",[interior_id])
  add_tag!(labels,"boundary",boundary_ids)
end

function _generate_pre_geolabel!(
  face_to_geolabel,
  face_to_cells,
  cell_to_faces,
  offset,
  interior_id,
  boundary_id,d,D)

  nfaces = length(face_to_cells)
  fill!(face_to_geolabel,interior_id)

  max_ncells_around = 2^(D-d)

  _generate_pre_geolabel_kernel!(
    face_to_geolabel,
    face_to_cells.data,
    face_to_cells.ptrs,
    cell_to_faces.data,
    cell_to_faces.ptrs,
    offset,boundary_id,max_ncells_around)

  face_to_geolabel
end

function _generate_pre_geolabel_kernel!(
  face_to_geolabel,
  face_to_cells_data,
  face_to_cells_ptrs,
  cell_to_faces_data,
  cell_to_faces_ptrs,
  offset,boundary_id,max_ncells_around)
  nfaces = length(face_to_geolabel)
  for face in 1:nfaces
    a = face_to_cells_ptrs[face]-1
    ncells_around = face_to_cells_ptrs[face+1] - (a+1)
    if ncells_around == 1
      icell_around = 1
      cell = face_to_cells_data[a+icell_around]
      b = cell_to_faces_ptrs[cell]-1
      nlfaces = cell_to_faces_ptrs[cell+1] - (b+1)
      for lface in 1:nlfaces
        face2 = cell_to_faces_data[b+lface]
        if face == face2
          face_to_geolabel[face] = lface + offset
          break
        end
      end
    elseif ncells_around != max_ncells_around
      face_to_geolabel[face] = boundary_id
    end
  end
end

function _fix_dface_geolabels!(
  dface_to_geolabel,
  jface_to_geolabel,
  dface_to_jfaces_data,
  dface_to_jfaces_ptrs,
  interior_id,boundary_id)

  ndfaces = length(dface_to_jfaces_ptrs)-1
  for dface in 1:ndfaces
    if dface_to_geolabel[dface] != boundary_id
      continue
    end
    a = dface_to_jfaces_ptrs[dface]
    b = dface_to_jfaces_ptrs[dface+1]-1
    for p in a:b
      jface = dface_to_jfaces_data[p]
      geolabel = jface_to_geolabel[jface]
      if geolabel != interior_id && geolabel != boundary_id
        dface_to_geolabel[dface] = geolabel
        break
      end
    end
  end
end

function _fill_subgrid_cartesian_face_labeling!(labels,topo,subdesc,desc,cmin)
   _fill_subgrid_cartesian_entities!(labels,topo,subdesc,desc,cmin)
   _add_cartesian_tags!(labels,topo)
end

function _fill_subgrid_cartesian_entities!(labels, topo, subdesc, desc, cmin)
  D = num_cell_dims(topo)
  d_to_dface_to_entity = labels.d_to_dface_to_entity
  gcis = CartesianIndices(desc.partition)
  subcis = CartesianIndices(subdesc.partition)

  polytope = first(get_polytopes(topo))
  face_labeling = labels
  offsets = get_offsets(polytope)
  interior_id = num_faces(polytope)
  boundary_id = -1

  minus_one_ci = CartesianIndex(tfill(-1, Val{D}()))

  polytope_d_face_to_jfaces = Matrix{Vector{Vector{Int}}}(undef, (D, D))
  for d = 0:(D-1)
    for j = d+1:D-1
      polytope_d_face_to_jfaces[d+1, j+1] = get_faces(polytope, d, j)
    end
  end

  face_deltas = _find_ncube_face_neighbor_deltas(polytope)
  for d = 0:(D-1)
    face_to_cells = get_faces(topo, d, D)
    cell_to_faces = get_faces(topo, D, d)
    face_to_geolabel = face_labeling.d_to_dface_to_entity[d+1]
    _generate_subgrid_pregeo_label!(
      face_to_geolabel,
      d,
      D,
      offsets,
      polytope_d_face_to_jfaces,
      minus_one_ci,
      face_to_cells,
      cell_to_faces,
      subcis,
      gcis,
      face_deltas,
      cmin,
      interior_id,
      boundary_id,
    )
  end
  _fix_geolabels(D, topo, d_to_dface_to_entity, interior_id, boundary_id)
  fill!(d_to_dface_to_entity[end], interior_id)
end

function _generate_subgrid_pregeo_label!(
  face_to_geolabel,
  d,
  D,
  offsets,
  polytope_d_face_to_jfaces,
  minus_one_ci,
  face_to_cells,
  cell_to_faces,
  subcis,
  gcis,
  face_deltas,
  cmin,
  interior_id,
  boundary_id,
)
  for face_gid = 1:length(face_to_geolabel)
    cell_gid = face_to_cells.data[face_to_cells.ptrs[face_gid]]
    a = cell_to_faces.ptrs[cell_gid]
    b = cell_to_faces.ptrs[cell_gid+1] - 1
    gci = (subcis[cell_gid] + minus_one_ci) + cmin

    face_lid = -1
    for j = a:b
      if (cell_to_faces.data[j] == face_gid)
        face_lid = j - a + 1
        break
      end
    end
    @assert face_lid != -1

    face_lid += offsets[d+1]
    # Check whether cell neighbour across face face_lid belongs to the
    # global grid. If yes, the current face is actually at the interior
    if ((gci + face_deltas[face_lid]) in gcis)
      face_to_geolabel[face_gid] = interior_id
    else
      cell_found = false
      for j = d+1:D-1
        dface_to_jfaces =
          polytope_d_face_to_jfaces[d+1, j+1][face_lid-offsets[d+1]]
        cell_found = _is_there_interior_cell_across_higher_dim_faces(
          dface_to_jfaces,
          offsets[j+1],
          gcis,
          gci,
          face_deltas,
        )
        cell_found && break
      end
      if (cell_found)
        # The current face is at least in two subgrid cells
        face_to_geolabel[face_gid] = boundary_id
      else
        # The current face is only in one subgrid cell
        face_to_geolabel[face_gid] = face_lid
      end
    end
  end
end


function _is_there_interior_cell_across_higher_dim_faces(
  dface_to_jfaces,
  offset_j,
  gcis,
  gci,
  face_deltas,
)
  for k in dface_to_jfaces
    jface_lid = k + offset_j
    if (isassigned(face_deltas, jface_lid))
      if ((gci + face_deltas[jface_lid]) in gcis)
        return true
      end
    end
  end
  return false
end

"""
  _find_ncube_face_neighbor_deltas(p::ExtrusionPolytope{D}) -> Vector{CartesianIndex}

  Given an n-cube type ExtrusionPolytope{D}, returns V=Vector{CartesianIndex} with as many
  entries as the number of faces in the boundary of the Polytope. For an entry face_lid
  in this vector, V[face_lid] returns what has to be added to the CartesianIndex of a
  cell in order to obtain the CartesianIndex of the cell neighbour of K across the face F
  with local ID face_lid.
"""
function _find_ncube_face_neighbor_deltas(p::ExtrusionPolytope{D}) where {D}
  nfaces = num_faces(p)
  delta_faces = Vector{CartesianIndex{D}}(undef, nfaces - 1)
  for face_lid = 1:nfaces-1
    delta_faces[face_lid] = _find_ncube_face_neighbor_delta(p, face_lid)
  end
  delta_faces
end

function _find_ncube_face_neighbor_delta(p::ExtrusionPolytope{D}, face_lid) where {D}
  @assert is_n_cube(p)
  result = fill(0, D)
  face = p.dface.nfaces[face_lid]
  for d = 1:D
    if (face.extrusion[d] == 0)
      result[d] = (face.anchor[d] == 0) ? -1 : 1
    end
  end
  return CartesianIndex(Tuple(result))
end

# Cartesian grid topology with periodic BC

function _cartesian_grid_topology_with_periodic_bcs(grid::UnstructuredGrid,
  isperiodic::NTuple,
  partition)

  cell_to_vertices, vertex_to_node =
    _generate_cell_to_vertices_from_grid(grid, isperiodic, partition)
  _generate_grid_topology_from_grid(grid,cell_to_vertices,vertex_to_node)
end

function _generate_cell_to_vertices_from_grid(grid::UnstructuredGrid,
  isperiodic::NTuple, partition)

  if is_first_order(grid)
    nodes = get_cell_nodes(grid)
    cell_to_vertices = copy(nodes)

    nnodes = num_nodes(grid)
    num_nodes_x_dir = [partition[i]+1 for i in 1:length(partition)]
    point_to_isperiodic, slave_point_to_point, slave_point_to_master_point =
      _generate_slave_to_master_point(num_nodes_x_dir,isperiodic, nnodes)

    vertex_to_point = findall( .! point_to_isperiodic)
    point_to_vertex = fill(-1,length(point_to_isperiodic))
    point_to_vertex[vertex_to_point] = 1:length(vertex_to_point)
    point_to_vertex[slave_point_to_point] = point_to_vertex[slave_point_to_master_point]

    cell_to_vertices = Table(lazy_map(Broadcasting(Reindex(point_to_vertex)),nodes))

    vertex_to_node = vertex_to_point
    node_to_vertex = point_to_vertex
  else
    @notimplemented
  end
  (cell_to_vertices,vertex_to_node, node_to_vertex)
end

function _generate_slave_to_master_point(num_nodes_x_dir::Vector{Int},
  isperiodic::NTuple, num_nodes::Int)

  point_to_isperiodic = fill(false,num_nodes)

  slave_ijk_bounds = Array{Any,1}(undef,length(isperiodic))
  master_ijk_bounds = Array{Any,1}(undef,length(isperiodic))

  linear_indices = LinearIndices(Tuple(num_nodes_x_dir))
  periodic_dirs = findall(x->x==true, isperiodic)
  for periodic_dir in periodic_dirs
    for i in 1:length(isperiodic)
      if i == periodic_dir
        slave_ijk_bounds[i] = num_nodes_x_dir[i]
      else
        slave_ijk_bounds[i] = 1:num_nodes_x_dir[i]
      end
    end
    slave_ijk_set = CartesianIndices(Tuple(slave_ijk_bounds))
    point_to_isperiodic[linear_indices[slave_ijk_set]] .= true
  end

  slave_point_to_point = findall( point_to_isperiodic)
  slave_point_to_master_point = Array{Int32,1}(undef,length(slave_point_to_point))

  cartesian_indices = CartesianIndices(Tuple(num_nodes_x_dir))
  for (i,point) in enumerate(slave_point_to_point)
    ijk = collect(cartesian_indices[point].I)
    for i in periodic_dirs
      if ijk[i] == num_nodes_x_dir[i]
        ijk[i] = 1
      end
    end

    master_point_ijk = CartesianIndex(Tuple(ijk))
    slave_point_to_master_point[i] = linear_indices[master_point_ijk]
  end

  point_to_isperiodic, slave_point_to_point, slave_point_to_master_point
end
