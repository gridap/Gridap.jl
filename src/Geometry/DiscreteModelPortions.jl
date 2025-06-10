
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
    restrict(model::DiscreteModel, cell_to_parent_cell::AbstractVector{<:Integer})
    restrict(model::DiscreteModel, parent_cell_to_mask::AbstractArray{Bool})
"""
@inline function restrict(model::DiscreteModel,args...;kwargs...)
  DiscreteModelPortion(model,args...;kwargs...)
end

@inline restrict(model::DiscreteModel, ::IdentityVector) = model

"""
"""
function DiscreteModelPortion(model::DiscreteModel, cell_to_parent_cell::AbstractVector{<:Integer})
  topo = get_grid_topology(model)
  labels = get_face_labeling(model)
  grid_p = restrict(get_grid(model),cell_to_parent_cell)
  topo_p, d_to_dface_to_parent_dface = restrict(topo,cell_to_parent_cell)
  labels_p = restrict(labels,d_to_dface_to_parent_dface)
  model_p = DiscreteModel(grid_p,topo_p,labels_p)
  DiscreteModelPortion(model_p,model,d_to_dface_to_parent_dface)
end

function DiscreteModelPortion(model::DiscreteModel, parent_cell_to_mask::AbstractArray{Bool})
  cell_to_parent_cell = findall(collect1d(parent_cell_to_mask))
  DiscreteModelPortion(model,cell_to_parent_cell)
end

function DiscreteModelPortion(model::DiscreteModel, parent_cell_to_mask::AbstractVector{Bool})
  cell_to_parent_cell = findall(parent_cell_to_mask)
  DiscreteModelPortion(model,cell_to_parent_cell)
end

function DiscreteModelPortion(model::DiscreteModel,grid_p::GridPortion)
  topo = get_grid_topology(model)
  labels = get_face_labeling(model)
  cell_to_parent_cell = grid_p.cell_to_parent_cell
  topo_p, d_to_dface_to_parent_dface = restrict(topo,cell_to_parent_cell)
  labels_p = restrict(labels,d_to_dface_to_parent_dface)
  model_p = DiscreteModel(grid_p,topo_p,labels_p)
  DiscreteModelPortion(model_p,model,d_to_dface_to_parent_dface)
end
