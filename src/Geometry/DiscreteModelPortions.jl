
"""
"""
struct DiscreteModelPortion{Dc,Dp} <: DiscreteModel{Dc,Dp}
  model::DiscreteModel{Dc,Dp}
  parent_model::DiscreteModel{Dc,Dp}
  d_to_dface_to_parent_dface::Vector{Vector{Int}}
end

get_grid(model::DiscreteModelPortion) = get_grid(model.model)

get_grid_topology(model::DiscreteModelPortion) = get_grid_topology(model.model)

get_face_labeling(model::DiscreteModelPortion) = get_face_labeling(model.model)

get_face_to_parent_face(model::DiscreteModelPortion,d::Integer) = model.d_to_dface_to_parent_dface[d+1]

get_cell_to_parent_cell(model::DiscreteModelPortion) = get_face_to_parent_face(model,num_cell_dims(model))

get_parent_model(model::DiscreteModelPortion) = model.parent_model

"""
"""
function DiscreteModelPortion(model::DiscreteModel, cell_to_parent_cell::AbstractVector{<:Integer})
  topo = get_grid_topology(model)
  labels = get_face_labeling(model)
  grid_p =  GridPortion(get_grid(model),cell_to_parent_cell)
  topo_p, d_to_dface_to_parent_dface = _grid_topology_portion(topo,cell_to_parent_cell)
  labels_p = _setup_labels_p(labels,d_to_dface_to_parent_dface)
  model_p = DiscreteModel(grid_p,topo_p,labels_p)
  DiscreteModelPortion(model_p,model,d_to_dface_to_parent_dface)
end

function DiscreteModelPortion(model::DiscreteModel, cell_to_is_in::AbstractArray{Bool})
  cell_to_parent_cell = findall(collect1d(cell_to_is_in))
  DiscreteModelPortion(model,cell_to_parent_cell)
end

function DiscreteModelPortion(model::DiscreteModel, cell_to_is_in::AbstractVector{Bool})
  cell_to_parent_cell = findall(cell_to_is_in)
  DiscreteModelPortion(model,cell_to_parent_cell)
end

function DiscreteModelPortion(model::DiscreteModel,grid_p::GridPortion)
  topo = get_grid_topology(model)
  labels = get_face_labeling(model)
  cell_to_parent_cell = grid_p.cell_to_parent_cell
  topo_p, d_to_dface_to_parent_dface = _grid_topology_portion(topo,cell_to_parent_cell)
  labels_p = _setup_labels_p(labels,d_to_dface_to_parent_dface)
  model_p = DiscreteModel(grid_p,topo_p,labels_p)
  DiscreteModelPortion(model_p,model,d_to_dface_to_parent_dface)
end

function _grid_topology_portion(topo,cell_to_parent_cell)
  D = num_cell_dims(topo)

  n_m_to_parent_nface_to_parent_mface = [
    Table(get_faces(topo,n,m)) for n in 0:D, m in 0:D ]

  d_to_dface_to_parent_dface = [
    _setup_dface_to_parent_dface(
      Val{d}(),
      Val{D}(),
      n_m_to_parent_nface_to_parent_mface[D+1,d+1],
      num_faces(topo,d),
      cell_to_parent_cell) for d in 0:D]

  d_to_dface_to_vertices  = [
    _setup_connectivities_d(
    n_m_to_parent_nface_to_parent_mface[d+1,0+1],
    d_to_dface_to_parent_dface[d+1],
    d_to_dface_to_parent_dface[0+1],
    num_faces(topo,0)) for d in 0:D]

  D = num_cell_dims(topo)
  vertex_coordinates = get_vertex_coordinates(topo)[d_to_dface_to_parent_dface[0+1]]
  cell_type = get_cell_type(topo)[d_to_dface_to_parent_dface[D+1]]
  polytopes = get_polytopes(topo)
  orientation_style = OrientationStyle(topo)

  topo_p = UnstructuredGridTopology(
    vertex_coordinates,
    d_to_dface_to_vertices,
    cell_type,
    polytopes,
    orientation_style)

  topo_p, d_to_dface_to_parent_dface
end

function _setup_dface_to_parent_dface(
  ::Val{Dc},
  ::Val{Dc},
  parent_cell_to_parent_dface::Table,
  num_parent_dfaces,
  cell_to_parent_cell) where Dc
  dface_to_parent_dface=Int[]
  parent_dface_touched = fill(false,num_parent_dfaces)
  for parent_cell in cell_to_parent_cell
    pini = parent_cell_to_parent_dface.ptrs[parent_cell]
    pend = parent_cell_to_parent_dface.ptrs[parent_cell+1]-1
    for p in pini:pend
      parent_dface = parent_cell_to_parent_dface.data[p]
      if (!parent_dface_touched[parent_dface])
         push!(dface_to_parent_dface,parent_dface)
         parent_dface_touched[parent_dface] = true
      end
    end
  end
  dface_to_parent_dface
end

function _setup_dface_to_parent_dface(
  ::Val{Df},
  ::Val{Dc},
  parent_cell_to_parent_dface::Table,
  num_parent_dfaces,
  cell_to_parent_cell) where {Df,Dc}
  parent_dface_touched = fill(false,num_parent_dfaces)
  for parent_cell in cell_to_parent_cell
    pini = parent_cell_to_parent_dface.ptrs[parent_cell]
    pend = parent_cell_to_parent_dface.ptrs[parent_cell+1]-1
    for p in pini:pend
      parent_dface = parent_cell_to_parent_dface.data[p]
      parent_dface_touched[parent_dface] = true
    end
  end
  dface_to_parent_dface = findall(parent_dface_touched)
end

function _setup_connectivities_d(
  parent_nface_parent_mface::Table,
  nface_to_parent_nface,
  mface_to_parent_mface,
  num_parent_mfaces)

  parent_mface_to_mface = fill(-1,num_parent_mfaces)
  parent_mface_to_mface[mface_to_parent_mface] .= 1:length(mface_to_parent_mface)

  nface_to_mface = parent_nface_parent_mface[nface_to_parent_nface]
  for p in 1:length(nface_to_mface.data)
    parent_mface = nface_to_mface.data[p]
    mface = parent_mface_to_mface[parent_mface]
    @check mface > 0
    nface_to_mface.data[p] = mface
  end
  nface_to_mface
end

function _setup_labels_p(labels,d_to_dface_to_parent_dface)

  D = length(d_to_dface_to_parent_dface)-1

  d_to_dface_to_entity = [
    zeros(
      eltype(labels.d_to_dface_to_entity[d+1]),
      length(d_to_dface_to_parent_dface[d+1]))
    for  d in 0:D]

  for d in 0:D
    face_to_label = d_to_dface_to_entity[d+1]
    oface_to_label = labels.d_to_dface_to_entity[d+1]
    face_to_oface = d_to_dface_to_parent_dface[d+1]
    _update_labels_dim!(face_to_label,oface_to_label,face_to_oface)
  end

  tag_to_entities = copy(labels.tag_to_entities)
  tag_to_name     = copy(labels.tag_to_name)
  FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)
end

function _update_labels_dim!(face_to_label,oface_to_label,face_to_oface)
  for face in 1:length(face_to_label)
    oface = face_to_oface[face]
    label = oface_to_label[oface]
    face_to_label[face] = label
  end
end
