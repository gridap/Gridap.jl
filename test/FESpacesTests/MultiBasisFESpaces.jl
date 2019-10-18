
function MultiShapeConformingFESpace(
  reffes,#::IndexCellValue{RefFE{D,T},N},
  trian::Triangulation{D,Z},
  graph::GridGraph,
  labels::FaceLabels,
  diri_tags::Vector{Int}) where {D,Z,N}
  args = _setup_multishape_conforming_fe_fields(reffes,trian,graph,labels,diri_tags,D)
  T = field_type(reffes[1])
  ConformingFESpace{D,Z,T}(args...)
end

function _setup_multishape_conforming_fe_fields(reffes,trian,graph,labels,diri_tags,D)
  reffe = reffes[1]
  dim_to_nface_eqclass, nfree, ndiri  = F._generate_dim_to_nface_to_dofs(
    reffe, graph, labels, diri_tags)
  cellvefs_dim = [connections(graph,D,i) for i in 0:D]
  offset = length.(dim_to_nface_eqclass)
  for i in 2:length(offset)
    offset[i] += offset[i-1]
  end
  offset = tuple(offset[1:(end-1)]...)
  cellvefs = local_append(offset, cellvefs_dim...)
  dofs_all = append(dim_to_nface_eqclass...)
  cell_eqclass = F.CellEqClass(cellvefs, dofs_all, reffe)
  values = broadcast(shfbasis,reffes)
  values = [shfbasis(reffe)]
  shb = broadcast(shfbasis,reffes)
  # pointer = [ i == ncells(trian) ? 2 : 1 for i in 1:ncells(trian)]
  # shb = CompressedCellValue(values,pointer)
  # shb = ConstantCellValue(shfbasis(reffe), ncells(trian))
  phi = CellGeomap(trian)
  basis = attachgeomap(shb,phi)
  return dim_to_nface_eqclass, cell_eqclass, nfree, ndiri, diri_tags,
    reffes, trian, graph, labels, basis
end
