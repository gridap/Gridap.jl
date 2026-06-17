# FE Space whose elements's polytope faces do not own all the DOFs of it's ref-FE

function FESpace(
  model::DiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{<:Loop}};
  conformity=nothing,
  trian = Triangulation(model),
  labels = get_face_labeling(model),
  dirichlet_tags=Int[],
  dirichlet_masks=nothing,
  constraint=nothing,
  vector_type=nothing,
  scale_dof=false,
  global_meshsize=nothing
)

  conf = Conformity(testitem(cell_reffe),conformity)
  cell_fe = FESpaces.CellFE(model, cell_reffe, conf; scale_dof, global_meshsize)
  fespace = FESpace(
    model, cell_fe; trian, labels, dirichlet_tags, dirichlet_masks, constraint, vector_type
  )

  #@assert fespace.

  reffe = testitem(cell_reffe)
  shapefuns = get_shapefuns(reffe)
  L = num_indep_components(value_type(shapefuns))

  ctopo = get_grid_topology(model)
  patch_vertex_map = LoopPatchVerticesMap(ctopo)
  cell_to_vertices = get_cell_vertices(ctopo)
  cell_patch_vertices = Table(lazy_map(patch_vertex_map, cell_to_vertices))

  for k in eachindex(cell_patch_vertices)
    k_dofs_ids = @view fespace.cell_dofs_ids[k]
    k_loop_ids = @view cell_patch_vertices[k]
    k_dof = 1
    for v in k_loop_ids
      for ll in 1:L
        k_dofs_ids[k_dof] = L*(v-1) + ll
        k_dof += 1
      end
    end
  end

  fespace
end


