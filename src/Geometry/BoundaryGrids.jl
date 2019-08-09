module BoundaryGrids

using Gridap
using Gridap.Helpers

export BoundaryGrid
import Gridap: Triangulation
import Gridap: BoundaryTriangulation
import Gridap: points
import Gridap: cells
import Gridap: celltypes
import Gridap: cellorders

struct BoundaryGrid{D,Z} <: Grid{D,Z}
  grid::Grid{D,Z}
  descriptor::BoundaryDescriptor
  
  function BoundaryGrid(
    grid::Grid{D,Z}, descriptor::BoundaryDescriptor) where {D,Z}
    @assert D == Z + 1
    new{D,Z}(grid,descriptor)
  end

end

for op in (:points,:cells,:celltypes,:cellorders)
  @eval begin
    $op(g::BoundaryGrid) = $op(g.grid)
  end
end

function Triangulation(grid::BoundaryGrid)
  trian = Triangulation(grid.grid)
  BoundaryTriangulation(trian,grid.descriptor)
end

function BoundaryGrid(model::DiscreteModel,tags,icell::Int=1)
  _tags = _setup_tags(model,tags)
  BoundaryGrid(model,_tags,icell)
end

function BoundaryGrid(model::DiscreteModel,tags::Vector{Int},icell::Int=1)
  d = celldim(model)
  labels = FaceLabels(model)
  oldfacet_to_label = labels_on_dim(labels,d-1)
  tag_to_labels = labels.tag_to_labels
  nfacets = length(oldfacet_to_label)
  oldfacet_to_mask = fill(false,nfacets)
  _setup_mask!(oldfacet_to_mask,oldfacet_to_label,tag_to_labels,tags)
  fgrid = Grid(model,d-1)
  facet_to_oldfacet = findall(oldfacet_to_mask)
  grid = GridPortion(fgrid,facet_to_oldfacet)
  graph = GridGraph(model)
  oldfacet_to_cells = connections(graph,d-1,d)
  cell_to_oldfacets = connections(graph,d,d-1)
  oldfacet_to_cell = get_local_item(oldfacet_to_cells,icell)
  oldfacet_to_lfacet = find_local_index(oldfacet_to_cell, cell_to_oldfacets)
  facet_to_cell = reindex(oldfacet_to_cell, facet_to_oldfacet)
  facet_to_lfacet = reindex(oldfacet_to_lfacet,facet_to_oldfacet)
  cellgrid = Grid(model,d)
  cell_to_extrussion = celltypes(cellgrid)
  cell_to_polytope = _cell_to_polytope(cell_to_extrussion)
  trian = Triangulation(cellgrid)
  phi = CellGeomap(trian)
  descriptor = BoundaryDescriptor(
    facet_to_cell,facet_to_lfacet,cell_to_polytope,phi)
  BoundaryGrid(grid,descriptor)
end

function BoundaryTriangulation(model::DiscreteModel,tags,icell::Int=1)
  grid = BoundaryGrid(model,tags,icell)
  Triangulation(grid)
end

function _setup_mask!(facet_to_mask,facet_to_label,tag_to_labels,tags)
  nfacets = length(facet_to_mask)
  for facet in 1:nfacets
    for tag in tags
      for label in tag_to_labels[tag]
        if facet_to_label[facet] == label
          facet_to_mask[facet] = true
        end
      end
    end
  end
end

function _cell_to_polytope(cell_to_extrussion)
  @notimplemented
end

function _cell_to_polytope(cell_to_extrussion::ConstantCellValue)
  extrussion = cell_to_extrussion.value
  l = cell_to_extrussion.length
  poly = Polytope(extrussion)
  ConstantCellValue(poly,l)
end

_setup_tags(model,tags) = tags

function _setup_tags(model,name::String)
  _setup_tags(model,[name,])
end

function _setup_tags(model,names::Vector{String})
  [ tag_from_name(model,s) for s in names ]
end

end # module
