
struct DiscreteModelMock <: DiscreteModel{2,2}
  grid::GridMock
  topo::GridTopologyMock
  labels::FaceLabeling
  function DiscreteModelMock()
    grid = GridMock()
    topo = GridTopologyMock()
    labels = _init_labeling(topo)
    new(grid,topo,labels)
  end
end

function _init_labeling(topo)
  d_to_num_dfaces = [ num_faces(topo,d) for d in 0:num_dims(topo)]
  labels = FaceLabeling(d_to_num_dfaces)
  get_face_entity(labels,0) .= get_isboundary_face(topo,0) .+ 1
  get_face_entity(labels,1) .= get_isboundary_face(topo,1) .+ 1
  get_face_entity(labels,2) .= get_isboundary_face(topo,2) .+ 1
  add_tag!(labels,"interior",[1,])
  add_tag!(labels,"boundary",[2,])
  add_tag_from_tags!(labels,"all",["interior","boundary"])
  labels
end

get_grid(model::DiscreteModelMock) = model.grid

get_grid_topology(model::DiscreteModelMock) = model.topo

get_face_labeling(model::DiscreteModelMock) = model.labels

