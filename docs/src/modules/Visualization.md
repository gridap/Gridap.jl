
```@meta
CurrentModule = Gridap.Visualization
```

# Gridap.Visualization

The `Gridap.Visualization` module provides tools to export Gridap grids (such as `Triangulation`, `DiscreteModel`) and objects defined on those grids (such as `Function` and `CellField`) into VTK format. These files can then be opened and visualized in external software like [ParaView](https://www.paraview.org/).

## High-level API

The primary entry point for visualization is `writevtk`.

```@docs
writevtk
```

For simulations that change over time, Gridap provides the PVD (ParaView Data) format to group multiple VTK files into a single time series.

```@docs
createpvd
savepvd
```

## Advanced & Low-level API

For more control over the visualization process, you can use `createvtk` to prepare a VTK object without writing it immediately, or work directly with `VisualizationData`.

```@docs
create_vtk_file
write_vtk_file
VisualizationData
visualization_data
VisualizationGrid
visualization_grid
```

