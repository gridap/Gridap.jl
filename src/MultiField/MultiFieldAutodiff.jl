
function CellData.SkeletonCellFieldPair(
  cf_plus::Union{MultiFieldFEFunction,MultiFieldCellField}, 
  cf_minus::Union{MultiFieldFEFunction,MultiFieldCellField}
)
  cfs = map(SkeletonCellFieldPair,cf_plus,cf_minus)
  MultiFieldCellField(cfs)
end
