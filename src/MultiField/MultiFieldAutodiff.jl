
# For skeletons

function FESpaces._get_cell_dof_values(f::MultiFieldFEFunction,trian::SkeletonTriangulation)
  uhs = f.single_fe_functions
  blockmask = [ is_change_possible(get_triangulation(uh),trian) for uh in uhs ]
  active_block_ids = findall(blockmask)
  active_block_data = Any[ FESpaces._get_cell_dof_values(uhs[i],trian) for i in active_block_ids ]
  nblocks = length(uhs)
  lazy_map(BlockMap(nblocks,active_block_ids),active_block_data...)
end

function CellData.SkeletonCellFieldPair(V::MultiFieldFESpace,cell_values::LazyArray{<:Fill{BlockMap{1}}})
  cfs = map(SkeletonCellFieldPair,V.spaces,cell_values.args)
  MultiFieldCellField(cfs)
end
