struct BoundaryDiscreteModel{Dc,Dp} <:DiscreteModel{Dc,Dp}
  model::DiscreteModel{Dc,Dp}
  bgmodel::DiscreteModel
  trian::BoundaryTriangulation
end

function BoundaryDiscreteModel(
  ::Type{<:Polytope{Df}},
  bgmodel::DiscreteModel,
  bgface_to_mask::AbstractArray) where Df

  @assert Df < num_cell_dims(bgmodel)

  bgmodel_Df = DiscreteModel(Polytope{Df},bgmodel)
  model = DiscreteModelPortion(bgmodel_Df,bgface_to_mask)
  trian = BoundaryTriangulation(bgmodel,bgface_to_mask)
  BoundaryDiscreteModel(model,bgmodel,trian)
end

get_grid(model::BoundaryDiscreteModel) = get_grid(model.model)
get_grid_topology(model::BoundaryDiscreteModel) = get_grid_topology(model.model)
get_face_labeling(model::BoundaryDiscreteModel) = get_face_labeling(model.model)
get_triangulation(a::BoundaryDiscreteModel) = a.trian
Triangulation(a::BoundaryDiscreteModel) = a.trian
