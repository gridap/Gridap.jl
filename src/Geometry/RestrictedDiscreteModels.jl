"""
"""
function DiscreteModel(parent_model::DiscreteModel,cell_to_parent_cell::AbstractVector{<:Integer})
  RestrictedDiscreteModel(parent_model,cell_to_parent_cell)
end

function DiscreteModel(parent_model::DiscreteModel,parent_cell_to_mask::AbstractVector{Bool})
  RestrictedDiscreteModel(parent_model,cell_to_mask)
end

function DiscreteModel(parent_model::DiscreteModel,labels::FaceLabeling,tags)
  parent_cell_to_mask = get_face_mask(labels,tags,num_cell_dims(parent_model))
  DiscreteModel(parent_model,parent_cell_to_mask)
end

function DiscreteModel(parent_model::DiscreteModel,tags)
  labels = get_face_labeling(parent_model)
  DiscreteModel(parent_model,labels,tags)
end

"""
"""
struct RestrictedDiscreteModel{Dc,Dp} <: DiscreteModel{Dc,Dc}
  model::DiscreteModelPortion{Dc,Dp}
end

function RestrictedDiscreteModel(parent_model::DiscreteModel, cell_to_parent_cell::Vector{Int})
  _model = DiscreteModelPortion(parent_model,cell_to_parent_cell)
  RestrictedDiscreteModel(_model)
end

function RestrictedDiscreteModel(parent_model::DiscreteModel, parent_cell_to_mask::Vector{Bool})
  _model = DiscreteModelPortion(parent_model,parent_cell_to_mask)
  RestrictedDiscreteModel(_model)
end

get_grid(model::RestrictedDiscreteModel) = get_grid(model.model)

get_grid_topology(model::RestrictedDiscreteModel) = get_grid_topology(model.model)

get_face_labeling(model::RestrictedDiscreteModel) = get_face_labeling(model.model)

get_face_to_parent_face(model::RestrictedDiscreteModel,d::Integer) = get_face_to_parent_face(model.model,d)

get_cell_to_parent_cell(model::RestrictedDiscreteModel) = get_cell_to_parent_cell(model.model)

get_parent_model(model::RestrictedDiscreteModel) = get_parent_model(model.model)

function Triangulation(model::RestrictedDiscreteModel)
  parent_trian = Triangulation(get_parent_model(model))
  cell_to_parent_cell = get_cell_to_parent_cell(model)
  TriangulationPortion(parent_trian,cell_to_parent_cell)
end

function get_triangulation(model::RestrictedDiscreteModel)
  Triangulation(model)
end

function Triangulation(::Type{ReferenceFE{d}},model::RestrictedDiscreteModel) where d
  @notimplemented
end

function Triangulation(model::RestrictedDiscreteModel,cell_to_parent_cell::AbstractVector{<:Integer})
  @notimplemented
end

function Triangulation(model::RestrictedDiscreteModel,cell_to_mask::AbstractVector{Bool})
  @notimplemented
end

function BoundaryTriangulation(model::RestrictedDiscreteModel,face_to_mask::Vector{Bool},icell_around::Integer)
  d = num_cell_dims(model)-1
  face_to_parent_face = get_face_to_parent_face(model,d)
  parent_model = get_parent_model(model)
  num_parent_faces = num_faces(parent_model,d)
  parent_face_to_mask = fill(false,num_parent_faces)
  parent_face_to_mask[face_to_parent_face] .= face_to_mask
  BoundaryTriangulation(parent_model,parent_face_to_mask)
end

function InterfaceTriangulation(model_in::RestrictedDiscreteModel,model_out::RestrictedDiscreteModel)
  cells_in = get_cell_to_parent_cell(model_in)
  cells_out = get_cell_to_parent_cell(model_out)
  model = get_parent_model(model_in)
  InterfaceTriangulation(model,cells_in,cells_out)
end

