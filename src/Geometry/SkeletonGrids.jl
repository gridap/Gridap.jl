module SkeletonGrids

using Gridap

export SkeletonGrid
import Gridap: points
import Gridap: cells
import Gridap: celltypes
import Gridap: cellorders
import Gridap: Triangulation
import Gridap: SkeletonTriangulation

struct SkeletonGrid{D,Z} <: Grid{D,Z}
  grid::Grid{D,Z}
  descriptor1::BoundaryDescriptor
  descriptor2::BoundaryDescriptor
  
  function SkeletonGrid(
    grid::Grid{D,Z},
    descriptor1::BoundaryDescriptor,
    descriptor2::BoundaryDescriptor) where {D,Z}
    @assert D == Z + 1
    new{D,Z}(grid,descriptor1,descriptor2)
  end

end

for op in (:points,:cells,:celltypes,:cellorders)
  @eval begin
    $op(g::SkeletonGrid) = $op(g.grid)
  end
end

function Triangulation(grid::SkeletonGrid)
  trian = Triangulation(grid.grid)
  SkeletonTriangulation(trian,grid.descriptor1,grid.descriptor2)
end

function SkeletonGrid(model::DiscreteModel, tags::Vector{Int})
  cell1 = 1
  cell2 = 2
  bgrid1 = BoundaryGrid(model,tags,cell1)
  bgrid2 = BoundaryGrid(model,tags,cell2)
  grid = bgrid1.grid
  descriptor1 = bgrid1.descriptor
  descriptor2 = bgrid2.descriptor
  SkeletonGrid(grid,descriptor1,descriptor2)
end

function SkeletonTriangulation(model::DiscreteModel,tags::Vector{Int})
  grid = SkeletonGrid(model,tags)
  Triangulation(grid)
end

end # module
