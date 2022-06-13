function MappedGrid(grid::Grid{Dc,Dp},uh::FEFunction) where {Dc,Dp}
  MappedGrid(grid,get_data(uh))
end

function MappedDiscreteModel(model::DiscreteModel{Dc,Dp},uh::FEFunction) where {Dc,Dp}
  MappedDiscreteModel(model,CellData.get_data(uh))
end
