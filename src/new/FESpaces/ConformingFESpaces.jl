
#  free_vals = zero_free_values(f)
#  dirichlet_vals = zero_dirichlet_values(f)
#  cache_vals = array_cache(cell_vals)
#  cache_dofs = array_cache(cell_dofs)
#
#  _free_and_dirichlet_values_fill!(
#    free_vals,
#    dirichlet_vals,
#    cache_vals,
#    cache_dofs,
#    cell_vals,
#    cell_dofs)
#
#function  _free_and_dirichlet_values_fill!(
#    free_vals,
#    dirichlet_vals,
#    cache_vals,
#    cache_dofs,
#    cell_vals,
#    cell_dofs)
#
#  for cell in 1:length(cell_vals)
#    vals = getindex!(cache_vals,cell_vals)
#    dofs = getindex!(cache_dofs,cell_dofs)
#    for (i,dof) in enumerate(dofs)
#      val = vals[i]
#      if dof > 0
#        free_vals[dof] = val
#      elseif dof < 0
#        dirichlet_vals[-dof] = val
#      else
#        @unreachable "dof ids either positibe or negative, not zero"
#      end
#    end
#  end
#
#end

