module BoundaryDescriptors

using Gridap

export BoundaryDescriptor

struct BoundaryDescriptor
  facet_to_cell::IndexCellNumber{<:Integer}
  facet_to_lfacet::IndexCellNumber{<:Integer}
  cell_to_polytope::IndexCellValue{<:Polytope}
  cell_phi::CellGeomap
  #TODO
  #facet_to_perm::IndexCellNumber{<:Integer}
end

end # module
