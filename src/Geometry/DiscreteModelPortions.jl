
"""
"""
struct DiscreteModelPortion{Dc,Dp} <: DiscreteModel{Dc,Dc}
  model::DiscreteModel{Dc,Dp}
  oldmodel::DiscreteModel{Dc,Dp}
  d_to_dface_to_old_dface::Vector{Vector{Int}}
end

get_grid(model::DiscreteModelPortion) = get_grid(model.model)

get_grid_topology(model::DiscreteModelPortion) = get_grid_topology(model.model)

get_face_labeling(model::DiscreteModelPortion) = get_face_labeling(model.model)

"""
"""
function DiscreteModelPortion(model::DiscreteModel, cell_to_oldcell::Vector{Int})
  grid_p =  GridPortion(get_grid(model),cell_to_oldcell)
  topo_p = GridTopology(grid_p)
  labels_p, d_to_dface_to_old_dface = _setup_labels_p(model,topo_p,cell_to_oldcell)
  model_p = DiscreteModel(grid_p,topo_p,labels_p)
  DiscreteModelPortion(model_p,model,d_to_dface_to_old_dface)
end

function _setup_labels_p(model,topo_p,cell_to_oldcell)

  topo = get_grid_topology(model)
  labels = get_face_labeling(model)
  D = num_cell_dims(model)

  T = eltype(cell_to_oldcell)
  d_to_dface_to_old_dface = Vector{Vector{T}}(undef,D+1)
  for d in 0:(D-1)
    cell_lface_to_face = Table(get_faces(topo_p,D,d))
    ocell_lface_to_oface = Table(get_faces(topo,D,d))
    nfaces = num_faces(topo_p,d)
    d_to_dface_to_old_dface[d+1] = _find_face_to_oldface(
      cell_lface_to_face,
      ocell_lface_to_oface,
      cell_to_oldcell,
      nfaces)
  end
  d_to_dface_to_old_dface[D+1] = cell_to_oldcell


  d_to_dface_to_entity = [
    zeros(eltype(labels.d_to_dface_to_entity[d+1]),num_faces(topo_p,d)) for  d in 0:D]
  for d in 0:D
    face_to_label = d_to_dface_to_entity[d+1]
    oface_to_label = labels.d_to_dface_to_entity[d+1]
    face_to_oface = d_to_dface_to_old_dface[d+1]
    _update_labels_dim!(face_to_label,oface_to_label,face_to_oface)
  end

  labels_p = FaceLabeling(d_to_dface_to_entity,labels.tag_to_entities,labels.tag_to_name)

  labels_p, d_to_dface_to_old_dface
end

function _find_face_to_oldface(
  cell_lface_to_face::Table, ocell_lface_to_oface::Table, cell_to_ocell, nfaces)

  # Assumes that the local faces have the same order in new and old

  T = eltype(cell_to_ocell)
  face_to_oface = zeros(T,nfaces)

  for cell in 1:length(cell_lface_to_face)
    ocell = cell_to_ocell[cell]
    a = cell_lface_to_face.ptrs[cell]-1
    b = cell_lface_to_face.ptrs[cell+1]
    nlfaces = b - (a+1)
    c = ocell_lface_to_oface.ptrs[ocell]-1
    for lface in 1:nlfaces
      face = cell_lface_to_face.data[a+lface]
      oface = ocell_lface_to_oface.data[c+lface]
      face_to_oface[face] = oface
    end
  end

  face_to_oface
end

function _update_labels_dim!(face_to_label,oface_to_label,face_to_oface)
  for face in 1:length(face_to_label)
    oface = face_to_oface[face]
    label = oface_to_label[oface]
    face_to_label[face] = label
  end
end

