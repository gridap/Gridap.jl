
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

function PatchFESpace(
  space::SingleFieldFESpace, ptopo::PatchTopology
)
  T = get_dof_value_type(space)
  model = get_background_model(get_triangulation(space))
  trian = PatchTriangulation(model,ptopo)
  fe_basis = GenericCellField(
    Fill(Fields.MockFieldArray(zeros(T,1)), num_cells(trian)), trian, ReferenceDomain()
  )
  return PatchFESpace(space, trian, ptopo, fe_basis)
end

function PatchFESpace(
  space::SingleFieldFESpace, trian::Triangulation, ptopo::PatchTopology, fe_basis
)
  model = get_background_model(trian)
  patch_dof_ids = generate_patch_dof_ids(space,ptopo)
  patch_trian = PatchTriangulation(model,ptopo)
  cell_dof_ids = _move_patch_coeffs(trian, patch_trian, patch_dof_ids)
  return PatchFESpace(space, trian, cell_dof_ids, fe_basis, ptopo)
end

function PatchFESpace(
  space::TrialFESpace, trian::Triangulation, 
  cell_dof_ids::AbstractArray, fe_basis,
  metadata = nothing
)
  test_space = PatchFESpace(space.space,trian,cell_dof_ids,fe_basis,metadata)
  TrialFESpace(space.dirichlet_values,test_space)
end

function PatchFESpace(
  space::SingleFieldFESpace, trian::Triangulation, 
  cell_dof_ids::AbstractArray, fe_basis,
  metadata = nothing
)
  cell_dof_ids = Table(cell_dof_ids)

  nfree = num_free_dofs(space)
  ndirichlet = num_dirichlet_dofs(space)
  dirichlet_dof_tag = get_dirichlet_dof_tag(space)
  ntags = num_dirichlet_tags(space)
  vector_type = get_vector_type(space)

  cell_is_dirichlet = collect(Bool,lazy_map(I -> any(i -> i < 0, I), cell_dof_ids))
  dirichlet_cells = collect(Int32,findall(cell_is_dirichlet))

  fe_dof_basis = CellData.MockCellDof(trian, DomainStyle(fe_basis))

  return UnconstrainedFESpace(
    vector_type, nfree, ndirichlet,
    cell_dof_ids, fe_basis, fe_dof_basis,
    cell_is_dirichlet, dirichlet_dof_tag, dirichlet_cells, ntags,
    metadata
  )
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
  space::SingleFieldFESpace, ptrian::PatchTriangulation,
)
  pface_to_dofs = get_cell_dof_ids(space,ptrian)
  patch_to_dofs = generate_patch_dof_ids(space,ptrian.ptopo)
  pface_to_patch = ptrian.glue.tface_to_patch
  pface_to_pdofs = find_local_index(pface_to_dofs,pface_to_patch,patch_to_dofs)
  return pface_to_pdofs
end

# Changes of domain

function CellData.change_domain(
  a::CellField,strian::PatchTriangulation,::ReferenceDomain,ttrian::PatchTriangulation,::ReferenceDomain
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
  a::CellField,strian::PatchTriangulation,::PhysicalDomain,ttrian::PatchTriangulation,::PhysicalDomain
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

function get_cell_fe_data(fun,f::SingleFieldFESpace,ttrian::PatchTriangulation)
  strian = get_triangulation(f)
  get_cell_fe_data(fun,f,strian,ttrian)
end

function get_cell_fe_data(fun,f,strian::Triangulation,ttrian::PatchTriangulation)
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

function get_cell_fe_data(fun,f,strian::PatchTriangulation,ttrian::PatchTriangulation)
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
  idegree::Integer = 2*order,
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
    order,idegree,metadata
  )
end
