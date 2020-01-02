
"""
    struct CartesianDiscreteModel{D,T,F} <: DiscreteModel{D,D}
      # Private Fields
    end
"""
struct CartesianDiscreteModel{D,T,F} <: DiscreteModel{D,D}
  grid::CartesianGrid{D,T,F}
  grid_topology::UnstructuredGridTopology{D,D,T,true}
  face_labeling::FaceLabeling
  @doc """
      CartesianDiscreteModel(desc::CartesianDescriptor)

  Inner constructor
  """
  function CartesianDiscreteModel(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
    grid = CartesianGrid(desc)
    _grid = UnstructuredGrid(grid)
    topo = UnstructuredGridTopology(_grid)
    nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
    labels = FaceLabeling(nfaces)
    _fill_cartesian_face_labeling!(labels,topo)
    new{D,T,F}(grid,topo,labels)
  end
end

"""
    CartesianDiscreteModel(args...)

Same args needed to construct a `CartesianDescriptor`
"""
function CartesianDiscreteModel(args...)
  desc = CartesianDescriptor(args...)
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
  face_nodes::Table{Int,Int32} = compute_face_nodes(model,d)
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
  _fill_cartesian_entitties!(labels,topo)
  _add_cartesian_tags!(labels,topo)
end

function _fill_cartesian_entitties!(labels,topo)
  D = num_cell_dims(topo)
  d_to_dface_to_entity = labels.d_to_dface_to_entity
  polytope = first(get_polytopes(topo))
  dim_to_offset = get_offsets(polytope)
  interior_id = num_faces(polytope)+1
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
  for d in 0:(D-2)
    for j in (d+1):(D-1)
      dface_to_jfaces = get_faces(topo,d,j)
      dface_to_geolabel = d_to_dface_to_entity[d+1]
      jface_to_geolabel = d_to_dface_to_entity[j+1]
      _fix_dface_geolabels!(
        dface_to_geolabel,
        jface_to_geolabel,
        dface_to_jfaces.data,
        dface_to_jfaces.ptrs,
        interior_id,boundary_id)
    end
  end
  fill!(d_to_dface_to_entity[end],interior_id)
end

function _add_cartesian_tags!(labels,topo)
  D = num_cell_dims(topo)
  polytope = first(get_polytopes(topo))
  interior_id = num_faces(polytope)+1
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
