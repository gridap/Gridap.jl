
# This FESpace preserves the order of the patch dofs between the original space 
# and this one. 
function FESpaceWithoutBCs(space::SingleFieldFESpace)
  nfree = num_free_dofs(space)
  ndir = num_dirichlet_dofs(space)
  free_ids = Int32(ndir+1):Int32(nfree+ndir)
  dir_ids = Int32(ndir):Int32(-1):Int32(1) # Reverse order
  return reindex_free_and_dirichlet_dof_ids(space,free_ids,dir_ids)
end

FESpaceWithoutBCs(space::ConstantFESpace) = space
FESpaceWithoutBCs(space::TrialFESpace) = FESpaceWithoutBCs(space.space)

function PatchFESpace(space::SingleFieldFESpace, ptopo::PatchTopology)
  vector_type = get_vector_type(space)

  nfree = num_free_dofs(space)
  cell_dof_ids = generate_patch_dof_ids(space,ptopo)

  model = get_background_model(get_triangulation(space))
  ptrian = Geometry.PatchTriangulation(model,ptopo)
  Dc = num_cell_dims(model)

  patch_to_ndofs = collect(lazy_map(length,cell_dof_ids))
  type_to_ndofs = unique(patch_to_ndofs)
  patch_to_type = collect(Int8,indexin(patch_to_ndofs,type_to_ndofs))
  type_to_basis = [MockFieldArray(zeros(Float64,ndofs)) for ndofs in type_to_ndofs]
  type_to_dofs = [ReferenceFEs.MockDofBasis(zeros(VectorValue{Dc,Float64},ndofs)) for ndofs in type_to_ndofs] 
  cell_shapefuns = expand_cell_data(type_to_basis,patch_to_type)
  cell_dof_basis = expand_cell_data(type_to_dofs,patch_to_type)
  fe_basis = GenericCellField(cell_shapefuns,ptrian,ReferenceDomain())
  fe_dof_basis = CellDof(cell_dof_basis,ptrian,ReferenceDomain())

  ndirichlet = num_dirichlet_dofs(space)
  dirichlet_dof_tag = get_dirichlet_dof_tag(space)
  ntags = num_dirichlet_tags(space)

  cell_is_dirichlet = collect(Bool,lazy_map(I -> any(i -> i < 0, I), cell_dof_ids))
  dirichlet_cells = collect(Int32,findall(cell_is_dirichlet))

  metadata = ptopo
  UnconstrainedFESpace(
    vector_type,
    nfree,
    ndirichlet,
    cell_dof_ids,
    fe_basis,
    fe_dof_basis,
    cell_is_dirichlet,
    dirichlet_dof_tag,
    dirichlet_cells,
    ntags,
    metadata
  )
end

function PatchFESpace(space::TrialFESpace,ptopo::PatchTopology)
  TrialFESpace(space.dirichlet_values,PatchFESpace(space.space,ptopo))
end

function generate_patch_dof_ids(
  space::SingleFieldFESpace, ptopo::PatchTopology
)
  trian = get_triangulation(space)
  Df = num_cell_dims(trian)
  glue = get_glue(trian,Val(Df))
  patch_to_lpface_to_mface = get_patch_faces(ptopo,Df)

  tface_to_dofs = get_cell_dof_ids(space)
  mface_to_dofs = extend(tface_to_dofs,glue.mface_to_tface)
  patch_to_dofs = Arrays.merge_entries(
    mface_to_dofs,patch_to_lpface_to_mface;
    acc = SortedSet{Int32}()
  )
  return patch_to_dofs
end

function generate_pface_to_pdofs(
  space::SingleFieldFESpace, ptopo::PatchTopology
)
  Df = num_cell_dims(get_triangulation(space))
  pface_to_dofs = get_cell_dof_ids(space)
  patch_to_dofs = generate_patch_dof_ids(space,ptopo)
  pface_to_patch = Geometry.get_pface_to_patch(ptopo,Df)
  pface_to_pdofs = find_local_index(pface_to_dofs,pface_to_patch,patch_to_dofs)
  return pface_to_pdofs
end

function generate_pface_to_pdofs(
  space::SingleFieldFESpace, ptrian::Geometry.PatchTriangulation,
)
  pface_to_dofs = get_cell_dof_ids(space,ptrian)
  patch_to_dofs = generate_patch_dof_ids(space,ptrian.ptopo)
  pface_to_patch = ptrian.glue.tface_to_patch
  pface_to_pdofs = find_local_index(pface_to_dofs,pface_to_patch,patch_to_dofs)
  return pface_to_pdofs
end

# Changes of domain

function CellData.change_domain(
  a::CellField,strian::Geometry.PatchTriangulation,::ReferenceDomain,ttrian::Geometry.PatchTriangulation,::ReferenceDomain
)
  if strian === ttrian
    return a
  end
  if is_change_possible(strian.trian,ttrian.trian)
    b = change_domain(a,strian.trian,ReferenceDomain(),ttrian.trian,ReferenceDomain())
    return CellData.similar_cell_field(b,CellData.get_data(b),ttrian,ReferenceDomain())
  end
  @assert num_cell_dims(strian) == num_cell_dims(ttrian)
  sglue = Geometry.get_patch_glue(strian)
  tglue = Geometry.get_patch_glue(ttrian)
  @check is_change_possible(sglue,tglue)
  return CellData.change_domain_ref_ref(a,ttrian,sglue,tglue)
end

function CellData.change_domain(
  a::CellField,strian::Geometry.PatchTriangulation,::PhysicalDomain,ttrian::Geometry.PatchTriangulation,::PhysicalDomain
)
  if strian === ttrian
    return a
  end
  if is_change_possible(strian.trian,ttrian.trian)
    b = change_domain(a,strian.trian,PhysicalDomain(),ttrian.trian,PhysicalDomain())
    return CellData.similar_cell_field(b,CellData.get_data(b),ttrian,PhysicalDomain())
  end
  @assert num_cell_dims(strian) == num_cell_dims(ttrian)
  sglue = Geometry.get_patch_glue(strian)
  tglue = Geometry.get_patch_glue(ttrian)
  @check is_change_possible(sglue,tglue)
  return CellData.change_domain_phys_phys(a,ttrian,sglue,tglue)
end

function get_cell_fe_data(fun,f::SingleFieldFESpace,ttrian::Geometry.PatchTriangulation)
  strian = get_triangulation(f)
  get_cell_fe_data(fun,f,strian,ttrian)
end

function get_cell_fe_data(fun,f,strian::Triangulation,ttrian::Geometry.PatchTriangulation)
  sface_to_data = fun(f)
  if strian === ttrian
    return sface_to_data
  end
  @check is_change_possible(strian,ttrian)
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  get_cell_fe_data(fun,sface_to_data,sglue,tglue)
end

function get_cell_fe_data(fun,f,strian::Geometry.PatchTriangulation,ttrian::Geometry.PatchTriangulation)
  sface_to_data = fun(f)
  if strian === ttrian
    return sface_to_data
  end
  if is_change_possible(strian.trian,ttrian.trian)
    D = num_cell_dims(strian)
    sglue = get_glue(strian,Val(D))
    tglue = get_glue(ttrian,Val(D))
    return get_cell_fe_data(fun,sface_to_data,sglue,tglue)
  end
  @assert num_cell_dims(strian) == num_cell_dims(ttrian)
  sglue = Geometry.get_patch_glue(strian)
  tglue = Geometry.get_patch_glue(ttrian)
  @check is_change_possible(sglue,tglue)
  return get_cell_fe_data(fun,sface_to_data,sglue,tglue)
end

############################################################################################

function PatchFESpace(model::DiscreteModel, ptopo::PatchTopology, args...; kwargs...)
  PatchFESpace(
    model, PatchTriangulation(model,ptopo), args...; kwargs...
  )
end

function PatchFESpace(trian::PatchTriangulation, args...; kwargs...)
  PatchFESpace(
    get_background_model(trian), trian, args...; kwargs...
  )
end

function PatchFESpace(
  model::DiscreteModel,trian::PatchTriangulation,::Type{T},order::Integer;
  space = :P, 
  vector_type = nothing,
  hierarchical = false, 
  orthonormal = false,
  local_kernel = nothing,
  kwargs...
) where T
  patch_grid = Geometry.bounding_box_grid(model, trian; δmin=0.1)
  patch_shapefuns, domain_style = get_polytopal_cell_shapefuns(
    patch_grid, T, order; space, hierarchical, orthonormal, local_kernel, domain_style = PhysicalDomain()
  )
  @check isa(domain_style, PhysicalDomain)
  vtype = ifelse(!isnothing(vector_type),vector_type,Vector{_dof_type(T)})
  return PatchFESpace(
    vtype, model, trian, patch_grid, patch_shapefuns, order, domain_style; kwargs...
  )
end

function PatchFESpace(
  vector_type::Type,
  model::DiscreteModel,
  trian::PatchTriangulation,
  patch_grid::Grid,
  patch_shapefuns::AbstractArray,
  order::Integer,
  domain_style::DomainStyle;
  labels = get_face_labeling(model),
  dirichlet_tags = Int[],
  dirichlet_masks = nothing,
)
  patch_to_faces = trian.glue.patch_to_tfaces
  face_to_patch = trian.glue.tface_to_patch

  npatches = num_cells(patch_grid)
  patch_to_ctype = Base.OneTo(npatches)
  ctype_to_ndofs = map(length, patch_shapefuns)

  ntags = length(dirichlet_tags)
  if ntags != 0
    face_to_tag = get_face_tag_index(labels,dirichlet_tags,Df)
    patch_to_tag = zeros(Int32, npatches)
    patch_is_dirichlet = zeros(Bool, npatches)
    for (patch,faces) in enumerate(patch_to_faces)
      tag = only(unique(face_to_tag[faces]))
      patch_to_tag[patch] = tag
      patch_is_dirichlet[patch] = !iszero(tag)
    end
    
    ctype_to_ldof_to_comp = lazy_map(n -> ones(Int16,n), ctype_to_ndofs)
    patch_dof_ids, nfree, ndir, dirichlet_dof_tag, dirichlet_patches = FESpaces.compute_discontinuous_cell_dofs(
      patch_to_ctype, ctype_to_ndofs, ctype_to_ldof_to_comp, patch_to_tag, dirichlet_masks
    )
  else
    ndir = 0
    dirichlet_dof_tag = Int8[]
    dirichlet_patches = Int32[]
    patch_is_dirichlet = fill(false,npatches)
    patch_dof_ids, nfree = FESpaces.compute_discontinuous_cell_dofs(patch_to_ctype,ctype_to_ndofs)
  end

  cell_shapefuns = lazy_map(Reindex(patch_shapefuns), face_to_patch)
  fe_basis = SingleFieldFEBasis(cell_shapefuns, trian, TestBasis(), domain_style)

  cell_dof_ids = Table(lazy_map(Broadcasting(Reindex(patch_dof_ids)), face_to_patch))
  cell_is_dirichlet = patch_is_dirichlet[face_to_patch]
  dirichlet_cells = reduce(vcat, (patch_to_faces[patch] for patch in dirichlet_patches); init = Int32[])

  metadata = nothing
  return PolytopalFESpace(
    vector_type,nfree,ndir,cell_dof_ids,fe_basis,
    cell_is_dirichlet,dirichlet_dof_tag,dirichlet_cells,ntags,
    order,metadata
  )
end
