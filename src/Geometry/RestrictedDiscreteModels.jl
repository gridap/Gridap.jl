
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

