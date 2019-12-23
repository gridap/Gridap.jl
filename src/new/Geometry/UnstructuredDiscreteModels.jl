
"""
    struct UnstructuredDiscreteModel{Dc,Dp,Tp,B} <: DiscreteModel{Dc,Dp}
      grid::UnstructuredGrid{Dc,Dp,Tp,B}
      grid_topology::UnstructuredGridTopology{Dc,Dp,Tp,B}
      face_labeling::FaceLabeling
    end
"""
struct UnstructuredDiscreteModel{Dc,Dp,Tp,B} <: DiscreteModel{Dc,Dp}
  grid::UnstructuredGrid{Dc,Dp,Tp,B}
  grid_topology::UnstructuredGridTopology{Dc,Dp,Tp,B}
  face_labeling::FaceLabeling
end

"""
    UnstructuredDiscreteModel(grid::Grid)
"""
function UnstructuredDiscreteModel(grid::Grid)
  _grid = UnstructuredGrid(grid)
  topo = UnstructuredGridTopology(_grid)
  nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
  labels = FaceLabeling(nfaces)
  UnstructuredDiscreteModel(_grid,topo,labels)
end

# Implementation of the interface

get_grid(model::UnstructuredDiscreteModel) = model.grid

get_grid_topology(model::UnstructuredDiscreteModel) = model.grid_topology

get_face_labeling(model::UnstructuredDiscreteModel) = model.face_labeling

