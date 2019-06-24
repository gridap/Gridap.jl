module BoundaryGrids

using Gridap

export BoundaryGrid
import Gridap: points
import Gridap: cells
import Gridap: celltypes
import Gridap: cellorders

struct BoundaryGrid{D,Z} <: Grid{D,Z}
  grid::Grid{D,Z}
  facet_to_cell::IndexCellNumber
  facet_to_lfacet::IndexCellNumber
  #TODO
  #facet_to_perm::IndexCellNumber
  
  function BoundaryGrid(
    grid::Grid{D,Z},
    facet_to_cell::IndexCellNumber,
    facet_to_lfacet::IndexCellNumber) where {D,Z}
    @assert D == Z + 1
    new{D,Z}(grid,facet_to_cell,facet_to_lfacet)
  end

end

for op in (:points,:cells,:celltypes,:cellorders)
  @eval begin
    $op(g::BoundaryGrid) = $op(g.grid)
  end
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
  graph = FullGridGraph(model) # TODO FullGridGraph not strictly needed
  oldfacet_to_cells = connections(graph,d-1,d)
  cell_to_oldfacets = connections(graph,d,d-1)
  oldfacet_to_cell = get_local_item(oldfacet_to_cells,icell)
  oldfacet_to_lfacet = find_local_index(oldfacet_to_cell, cell_to_oldfacets)
  facet_to_cell = reindex(oldfacet_to_cell, facet_to_oldfacet)
  facet_to_lfacet = reindex(oldfacet_to_lfacet,facet_to_oldfacet)
  BoundaryGrid(grid,facet_to_cell,facet_to_lfacet)
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

end # module
