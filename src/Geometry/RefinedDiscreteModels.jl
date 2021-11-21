# setp 1
function newest_vertex_bisection(grid::Grid,cell_mask::AbstractVector{<:Bool})
    #get_faces(top
    cell_coords = get_cell_coordinates(grid)
    @show typeof(cell_coords)
    #@show new_val = VectorValue{2, Float64}(1.5, 1.5, 1.5)
    @show new_val = [VectorValue(1.5, 1.5, 1.5)]
    cell_coords = lazy_append(cell_coords, new_val)
    @show cell_coords 
  # tod
  #ref_grid
end

# step 2
function newest_vertex_bisection(model::DiscreteModel,cell_mask::AbstractVector{<:Bool})
  grid  = get_grid(model)
  ref_grid = newest_vertex_bisection(grid,cell_mask)
  #ref_topo = GridTopology(grid)
  #labels = get_face_labelling(model)
  #ref_labels = # Compute them from the original labels (This is perhaps the most tedious part)
  #ref_model = DiscreteModel(ref_grid,ref_topo,ref_labels)
  #ref_model
end
