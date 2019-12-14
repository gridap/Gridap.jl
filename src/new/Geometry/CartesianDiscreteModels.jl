
"""
    struct CartesianDiscreteModel{D,T,F} <: DiscreteModel{D,D}
      # Private Fields
    end
"""
struct CartesianDiscreteModel{D,T,F} <: DiscreteModel{D,D}
  grid::CartesianGrid{D,T,F}
  model::UnstructuredDiscreteModel{D,D}
  @doc """
      CartesianDiscreteModel(desc::CartesianDescriptor)
  """
  function CartesianDiscreteModel(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
    grid = CartesianGrid(desc)
    model = UnstructuredDiscreteModel(grid)
    _fill_cartesian_face_labeling!(model)
    new{D,T,F}(grid,model)
  end
end

# non-default traits

OrientationStyle(::Type{<:CartesianDiscreteModel}) = Val{true}()

"""
    CartesianDiscreteModel(args...)

Same args needed to construct a `CartesianDescriptor`
"""
function CartesianDiscreteModel(args...)
  desc = CartesianDescriptor(args...)
  CartesianDiscreteModel(desc)
end

# Delegated to the grid
"""
    get_cartesian_descriptor(grid::CartesianDiscreteModel)

Get the descriptor of the Cartesian model
"""
get_cartesian_descriptor(a::CartesianDiscreteModel) = get_cartesian_descriptor(a.grid)

get_cell_type(g::CartesianDiscreteModel) = get_cell_type(g.grid)

get_cell_map(g::CartesianDiscreteModel) = get_cell_map(g.grid)

num_nodes(g::CartesianDiscreteModel) = num_nodes(g.grid)

# Delegated to the model

get_node_coordinates(g::CartesianDiscreteModel) = get_node_coordinates(g.model)

num_faces(g::CartesianDiscreteModel,d::Integer) = num_faces(g.model,d)

get_faces(g::CartesianDiscreteModel,dimfrom::Integer,dimto::Integer) = get_faces(g.model,dimfrom,dimto)

get_vertex_node(g::CartesianDiscreteModel) = get_vertex_node(g.model)

get_node_face_owner(g::CartesianDiscreteModel) = get_node_face_owner(g.model)

get_face_nodes(g::CartesianDiscreteModel,d::Integer) = get_face_nodes(g.model,d)

get_face_own_nodes(g::CartesianDiscreteModel,d::Integer) = get_face_own_nodes(g.model,d)

get_cell_nodes(g::CartesianDiscreteModel) = get_cell_nodes(g.model)

get_isboundary_face(g::CartesianDiscreteModel,d::Integer) = get_isboundary_face(g.model,d)

get_face_reffe_type(g::CartesianDiscreteModel,d::Integer) = get_face_reffe_type(g.model,d)

get_face_polytope_type(g::CartesianDiscreteModel,d::Integer) = get_face_polytope_type(g.model,d)

get_reffes(
 ::Type{<:ReferenceFE{d}},g::CartesianDiscreteModel) where d = get_reffes(NodalReferenceFE{d},g.model)

get_polytopes(
 ::Type{<:Polytope{d}},g::CartesianDiscreteModel) where d = get_polytopes(Polytope{d},g.model)

get_face_labeling(g::CartesianDiscreteModel) = get_face_labeling(g.model)

"""
    get_polytope(::Type{<:Polytope{d}},model::CartesianDiscreteModel) where d
"""
function get_polytope(::Type{<:Polytope{d}},model::CartesianDiscreteModel) where d
  polytopes = get_polytopes(Polytope{d},model)
  @assert length(polytopes) == 1
  first(polytopes)
end

"""
    get_polytope(model::CartesianDiscreteModel)
"""
function get_polytope(model::CartesianDiscreteModel)
  D = num_cell_dims(model)
  get_polytope(Polytope{D},model)
end


# Helpers

function _fill_cartesian_face_labeling!(model)
  _fill_cartesian_entitties!(model)
  _add_cartesian_tags!(model)
end

function _fill_cartesian_entitties!(model)
  D = num_cell_dims(model)
  labels = get_face_labeling(model)
  d_to_dface_to_entity = labels.d_to_dface_to_entity
  polytope = first(get_polytopes(Polytope{D},model))
  dim_to_offset = get_offsets(polytope)
  interior_id = num_faces(polytope)+1
  boundary_id = -1
  for d in 0:(D-1)
    face_to_cells = get_faces(model,d,D)
    cell_to_faces = get_faces(model,D,d)
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
      dface_to_jfaces = get_faces(model,d,j)
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

function _add_cartesian_tags!(model)
  D = num_cell_dims(model)
  labels = get_face_labeling(model)
  polytope = first(get_polytopes(Polytope{D},model))
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
