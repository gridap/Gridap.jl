
struct LinearizedFESpace{M,RM,LF} <: SingleFieldFESpace
  model::M
  refined_model::RM
  linear_fe_space::LF
end

function _linearize_reffe(reffe::Tuple{<:Lagrangian,Any,Any})
   (reffe[1],(reffe[2][1],1),reffe[3])
end

function _num_cells_x_dim(reffe::Tuple{<:Lagrangian,Any,Any})
   reffe[2][2]
end

function LinearizedFESpace(model::CartesianDiscreteModel{Dc},
                           reffe::Tuple{<:Lagrangian,Any,Any}; kwargs...) where Dc
  reffe_linearized=_linearize_reffe(reffe)
  ncells_x_dim=_num_cells_x_dim(reffe)
  ref_model=Adaptivity.refine(model,ncells_x_dim)
  linear_fe_space=FESpace(ref_model,reffe_linearized;kwargs...)
  LinearizedFESpace(model,ref_model,linear_fe_space)
end

# FESpace interface
FESpaces.get_free_dof_ids(f::LinearizedFESpace) = get_free_dof_ids(f.linear_fe_space)
FESpaces.zero_free_values(f::LinearizedFESpace) = zero_free_values(f.linear_fe_space)

function FESpaces.get_cell_dof_ids(f::LinearizedFESpace)
  get_cell_dof_ids(f.linear_fe_space)
end

function FESpaces.get_fe_dof_basis(f::LinearizedFESpace)
  get_fe_dof_basis(f.linear_fe_space)
end

Geometry.get_triangulation(f::LinearizedFESpace) = Triangulation(f.refined_model)
FESpaces.get_dof_value_type(f::LinearizedFESpace) = get_dof_value_type(f.linear_fe_space)
FESpaces.get_vector_type(f::LinearizedFESpace) = get_vector_type(f.linear_fe_space)
FESpaces.get_cell_is_dirichlet(f::LinearizedFESpace) = get_cell_is_dirichlet(f.linear_fe_space)
FESpaces.ConstraintStyle(::Type{<:LinearizedFESpace}) = FESpaces.UnConstrained()

# SingleFieldFESpace interface
FESpaces.get_dirichlet_dof_ids(f::LinearizedFESpace) = get_dirichlet_dof_ids(f.linear_fe_space)
FESpaces.num_dirichlet_tags(f::LinearizedFESpace) = num_dirichlet_tags(f.linear_fe_space)
FESpaces.zero_dirichlet_values(f::LinearizedFESpace) = zero_dirichlet_values(f.linear_fe_space)
FESpaces.get_dirichlet_dof_tag(f::LinearizedFESpace) = get_dirichlet_dof_tag(f.linear_fe_space)

function FESpaces.scatter_free_and_dirichlet_values(f::LinearizedFESpace,free_values,dirichlet_values)
  scatter_free_and_dirichlet_values(f.linear_fe_space,free_values,dirichlet_values)
end

function FESpaces.gather_free_and_dirichlet_values!(free_vals,dirichlet_vals,f::LinearizedFESpace,cell_vals)
  gather_free_and_dirichlet_values!(free_vals,dirichlet_vals,f.linear_fe_space,cell_vals)
end

function gather_dirichlet_values!(dirichlet_vals,f::LinearizedFESpace,cell_vals)
  gather_dirichlet_values!(dirichlet_vals,f.linear_fe_space,cell_vals)
end

function FESpaces.get_fe_basis(f::LinearizedFESpace)
  fe_basis = get_fe_basis(f.linear_fe_space)
end

# The following defitions are required to support the correct execution
# of the collect_cell_matrix(trial::FESpace,test::FESpace,a::DomainContribution)
# function for different combinations of regular and adapted triangulations
# corresponding to FESpace's triangulation and domain of integration
# Please, use with care, they are only valid for the most obvious cases.
function Gridap.FESpaces.get_cell_fe_data(fun,f,ttrian::AdaptedTriangulation)
  sface_to_data = fun(f)
  strian = get_triangulation(f)
  _get_cell_fe_data(sface_to_data, strian, ttrian)
end

# strian: triangulation of the FESpace, cell-wise related data (e.g., cell_dof_ids)
# ttrian: triangulation of the Measure
function _get_cell_fe_data(sface_to_data, strian::Triangulation, ttrian::Triangulation)
  _get_cell_fe_data_trian_trian_body(sface_to_data, strian, ttrian)
end

function _get_cell_fe_data_trian_trian_body(sface_to_data, strian, ttrian)
  if strian === ttrian
    return sface_to_data
  end
  @assert is_change_possible(strian,ttrian)
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  Gridap.FESpaces.get_cell_fe_data(sface_to_data,sglue,tglue)
end

function _get_cell_fe_data(sface_to_data, strian::AdaptedTriangulation, ttrian::Triangulation)
  _get_cell_fe_data_trian_trian_body(sface_to_data, strian, ttrian)
 end

# strian coarse
# ttrian refined
function _get_cell_fe_data(sface_to_data, strian::Triangulation, ttrian::AdaptedTriangulation)
  Gridap.Helpers.@check get_background_model(strian) === get_parent(get_adapted_model(ttrian))
  Gridap.Adaptivity.c2f_reindex(sface_to_data,get_adapted_model(ttrian).glue)
end

function _get_cell_fe_data(sface_to_data, strian::AdaptedTriangulation, ttrian::AdaptedTriangulation)
  Gridap.Helpers.@check get_background_model(strian)===get_background_model(ttrian)
  _get_cell_fe_data_trian_trian_body(sface_to_data, strian, ttrian)
end


function Gridap.FESpaces._compute_cell_ids(uh,ttrian)
  strian = get_triangulation(uh)
  Gridap.FESpaces._compute_cell_ids(uh, strian, ttrian)
end

function Gridap.FESpaces._compute_cell_ids(uh, strian::Triangulation, ttrian::Triangulation)
  _compute_cell_ids_body(uh,strian,ttrian)
end

function _compute_cell_ids_body(uh, strian, ttrian)
  if strian === ttrian
    return collect(IdentityVector(Int32(num_cells(strian))))
  end
  @check is_change_possible(strian,ttrian)
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  @notimplementedif !isa(sglue,FaceToFaceGlue)
  @notimplementedif !isa(tglue,FaceToFaceGlue)
  scells = IdentityVector(Int32(num_cells(strian)))
  mcells = extend(scells,sglue.mface_to_tface)
  tcells = lazy_map(Reindex(mcells),tglue.tface_to_mface)
  collect(tcells)
end

function Gridap.FESpaces._compute_cell_ids(uh, strian::Triangulation, ttrian::AdaptedTriangulation)
  Gridap.Helpers.@check get_background_model(strian) === get_parent(get_adapted_model(ttrian))
  return collect(IdentityVector(Int32(num_cells(ttrian))))
end

function Gridap.FESpaces._compute_cell_ids(uh, strian::AdaptedTriangulation, ttrian::Triangulation)
  Gridap.Helpers.@notimplemented
end

function Gridap.FESpaces._compute_cell_ids(uh,
                                           strian::AdaptedTriangulation,
                                           ttrian::AdaptedTriangulation)
  _compute_cell_ids_body(uh, strian, ttrian)
end

function Gridap.FESpaces._jacobian(f,uh,fuh::DomainContribution)
  terms = DomainContribution()
  for trian in get_domains(fuh)
    g = Gridap.FESpaces._change_argument(jacobian,f,trian,uh)
    cell_u = get_cell_dof_values(uh,trian)
    cell_id = Gridap.FESpaces._compute_cell_ids(uh,trian)
    cell_grad = Gridap.FESpaces.autodiff_array_jacobian(g,cell_u,cell_id)
    Gridap.CellData.add_contribution!(terms,trian,cell_grad)
  end
  terms
end
