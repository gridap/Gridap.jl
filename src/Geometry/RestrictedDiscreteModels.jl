"""
"""
function DiscreteModel(model::DiscreteModel,cell_to_oldcell::AbstractVector{<:Integer})
  RestrictedDiscreteModel(model,cell_to_oldcell)
end

function DiscreteModel(model::DiscreteModel,cell_to_mask::AbstractVector{Bool})
  RestrictedDiscreteModel(model,cell_to_mask)
end

function DiscreteModel(model::DiscreteModel,labels::FaceLabeling,tags)
  cell_to_mask = get_face_mask(labels,tags,num_cell_dims(model))
  DiscreteModel(model,cell_to_mask)
end

function DiscreteModel(model::DiscreteModel,tags)
  labels = get_face_labeling(model)
  DiscreteModel(model,labels,tags)
end




"""
"""
struct RestrictedDiscreteModel{Dc,Dp} <: DiscreteModel{Dc,Dc}
  model::DiscreteModelPortion{Dc,Dp}
end

function RestrictedDiscreteModel(model::DiscreteModel, cell_to_oldcell::Vector{Int})
  _model = DiscreteModelPortion(model,cell_to_oldcell)
  RestrictedDiscreteModel(_model)
end

function RestrictedDiscreteModel(model::DiscreteModel, cell_to_mask::Vector{Bool})
  _model = DiscreteModelPortion(model,cell_to_mask)
  RestrictedDiscreteModel(_model)
end

get_grid(model::RestrictedDiscreteModel) = get_grid(model.model)

get_grid_topology(model::RestrictedDiscreteModel) = get_grid_topology(model.model)

get_face_labeling(model::RestrictedDiscreteModel) = get_face_labeling(model.model)

get_face_to_oldface(model::RestrictedDiscreteModel,d::Integer) = get_face_to_oldface(model.model,d)

get_cell_to_oldcell(model::RestrictedDiscreteModel) = get_cell_to_oldcell(model.model)

get_oldmodel(model::RestrictedDiscreteModel) = get_oldmodel(model.model)

function RestrictedTriangulation(model::RestrictedDiscreteModel)
  oldtrian = Triangulation(get_oldmodel(model))
  cell_to_oldcell = get_cell_to_oldcell(model)
  RestrictedTriangulation(oldtrian,cell_to_oldcell)
end

function Triangulation(model::RestrictedDiscreteModel)
  RestrictedTriangulation(model)
end

function get_triangulation(model::RestrictedDiscreteModel)
  RestrictedTriangulation(model)
end

function Triangulation(::Type{ReferenceFE{d}},model::RestrictedDiscreteModel) where d
  @notimplemented
end

function Triangulation(model::RestrictedDiscreteModel,cell_to_oldcell::AbstractVector{<:Integer})
  @notimplemented
end

function Triangulation(model::RestrictedDiscreteModel,cell_to_mask::AbstractVector{Bool})
  @notimplemented
end

function BoundaryTriangulation(model::RestrictedDiscreteModel,face_to_mask::Vector{Bool},icell_around::Integer)
  d = num_cell_dims(model)-1
  face_to_oldface = get_face_to_oldface(model,d)
  oldmodel = get_oldmodel(model)
  num_oldfaces = num_faces(oldmodel,d)
  oldface_to_mask = fill(false,num_oldfaces)
  oldface_to_mask[face_to_oldface] .= face_to_mask
  BoundaryTriangulation(oldmodel,oldface_to_mask)
end

function InterfaceTriangulation(model_in::RestrictedDiscreteModel,model_out::RestrictedDiscreteModel)
  cells_in = get_cell_to_oldcell(model_in)
  cells_out = get_cell_to_oldcell(model_out)
  model = get_oldmodel(model_in)
  InterfaceTriangulation(model,cells_in,cells_out)
end

