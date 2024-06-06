
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
  order = reffe[2][2]
  ref_model=Adaptivity.refine(model,order)
  linear_fe_space=FESpace(ref_model,reffe_linearized;kwargs...)
  LinearizedFESpace(model,ref_model,linear_fe_space)
end

function LinearizedFESpace(model::CartesianDiscreteModel{Dc},
                           reffe::Tuple{<:Nedelec,Any,Any}; kwargs...) where Dc
  reffe_linearized=_linearize_reffe(reffe)
  order = reffe[2][2]
  ref_model=Adaptivity.refine(model,order-1)
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
  _get_cell_fe_data(fun,sface_to_data, strian, ttrian)
end

# strian: triangulation of the FESpace, cell-wise related data (e.g., cell_dof_ids)
# ttrian: triangulation of the Measure
function _get_cell_fe_data(fun, sface_to_data, strian::Triangulation, ttrian::Triangulation)
  _get_cell_fe_data_trian_trian_body(fun, sface_to_data, strian, ttrian)
end

function _get_cell_fe_data_trian_trian_body(fun, sface_to_data, strian, ttrian)
  if strian === ttrian
    return sface_to_data
  end
  @assert is_change_possible(strian,ttrian)
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  Gridap.FESpaces.get_cell_fe_data(fun,sface_to_data,sglue,tglue)
end

function _get_cell_fe_data(fun,sface_to_data, strian::AdaptedTriangulation, ttrian::Triangulation)
  _get_cell_fe_data_trian_trian_body(fun,sface_to_data, strian, ttrian)
 end

# strian coarse
# ttrian refined
function _get_cell_fe_data(fun, sface_to_data, strian::Triangulation, ttrian::AdaptedTriangulation)
  Gridap.Helpers.@check get_background_model(strian) === get_parent(get_adapted_model(ttrian))
  Gridap.Helpers.@check isa(get_adapted_model(ttrian).glue.n2o_faces_map[end],Vector)
  Dc = num_cell_dims(strian)
  adapt_glue=get_adapted_model(ttrian).glue
  sglue=get_glue(strian,Val(Dc))
  s_mface_to_tface=sglue.mface_to_tface
  tglue=get_glue(ttrian,Val(Dc))
  sface_to_data_reindexed=_reindex_sface_to_data(fun,sface_to_data,s_mface_to_tface,adapt_glue,tglue)
  _get_cell_fe_data_trian_trian_body(fun,sface_to_data_reindexed, ttrian.trian, ttrian)
end

function _reindex_sface_to_data(fun,sface_to_data,s_mface_to_tface,adapt_glue,tglue)
  t_tface_to_mface=tglue.tface_to_mface
  s_coarse_mfaces=adapt_glue.n2o_faces_map[end][t_tface_to_mface]
  lazy_map(Reindex(sface_to_data),s_mface_to_tface[s_coarse_mfaces])
end 

function _reindex_sface_to_data(fun,sface_to_data,s_mface_to_tface,adapt_glue,tglue::SkeletonPair)
  t_tface_to_mface_plus  = tglue.plus.tface_to_mface
  t_tface_to_mface_minus = tglue.minus.tface_to_mface
  s_coarse_mfaces_plus  = adapt_glue.n2o_faces_map[end][t_tface_to_mface_plus]
  s_coarse_mfaces_minus = adapt_glue.n2o_faces_map[end][t_tface_to_mface_minus]
  sface_to_data_plus_reindex = lazy_map(Reindex(sface_to_data),s_mface_to_tface[s_coarse_mfaces_plus])
  sface_to_data_minus_reindex = lazy_map(Reindex(sface_to_data),s_mface_to_tface[s_coarse_mfaces_minus])
  _combine_plus_minus(fun,sface_to_data_plus_reindex,sface_to_data_minus_reindex)
end 

function _combine_plus_minus(fun,p,m)
  lazy_map(BlockMap(2,[1,2]),p,m)
end 

function _combine_plus_minus(fun::typeof(get_cell_is_dirichlet),p,m)
  lazy_map((p,m)->(p||m),p,m)
end 

function _get_cell_fe_data(fun, 
                           sface_to_data, 
                           strian::AdaptedTriangulation, 
                           ttrian::AdaptedTriangulation)
  if (get_background_model(strian) === get_parent(get_adapted_model(ttrian)))
    sface_to_data_reindexed=Gridap.Adaptivity.o2n_reindex(sface_to_data,get_adapted_model(ttrian).glue)
    # Strictly speaking, the next line is not guaranteed to work in the most general case. 
    # I am assuming that Triangulation(get_background_model(ttrian)) is the proper one that 
    # matches sface_to_data_reindexed, but this is not true in general. For the particular (and most frequent)
    # case in which sface_to_data holds cell-wise data coming from a FESpace, this is true.
    _get_cell_fe_data_trian_trian_body(fun,sface_to_data_reindexed, Triangulation(get_background_model(ttrian)), ttrian)
  elseif (get_background_model(strian) === get_background_model(ttrian))
    _get_cell_fe_data_trian_trian_body(fun,sface_to_data, strian, ttrian)
  else
    Gridap.Helpers.@notimplemented
  end 
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

### BEGIN OPTIMIZATIONS
function Gridap.Adaptivity.get_n2o_reference_coordinate_map(
            g::AdaptivityGlue{Gridap.Adaptivity.RefinementGlue,Dc,A,B,<:FillArrays.Fill}) where {Dc,A,B}
  rrules    = Gridap.Adaptivity.get_new_cell_refinement_rules(g)
  cell_maps = Gridap.Geometry.get_cell_map(rrules.value)
  return lazy_map(Reindex(cell_maps),g.n2o_cell_to_child_id)
end

function Gridap.Adaptivity.get_n2o_reference_coordinate_map(
            g::AdaptivityGlue{Gridap.Adaptivity.RefinementGlue,Dc,A,B,<:CompressedArray}) where {Dc,A,B}
  
  rrules = Gridap.Adaptivity.get_new_cell_refinement_rules(g)
  if (length(rrules.values)==1)
     cell_maps = Gridap.Geometry.get_cell_map(rrules.values[1])
     return lazy_map(Reindex(cell_maps),g.n2o_cell_to_child_id)
  else
    @assert length(rrules.values)==2
    rrules_types=map(typeof,map(Gridap.Adaptivity.RefinementRuleType,rrules.values))
    @assert Gridap.Adaptivity.WithoutRefinement in rrules_types
    @assert Gridap.Adaptivity.GenericRefinement in rrules_types
    cell_maps_1=Gridap.Geometry.get_cell_map(rrules.values[1])
    cell_maps_2=Gridap.Geometry.get_cell_map(rrules.values[2])
    if (rrules_types[1]==Gridap.Adaptivity.WithoutRefinement)
      cell_maps=[cell_maps_2...,cell_maps_1...]
    else
      cell_maps=[cell_maps_1...,cell_maps_2...]
    end
    wo_refinement_map_id=maximum(g.n2o_cell_to_child_id)+1
    cell_maps_ptrs=Vector{Int}(undef,length(g.n2o_cell_to_child_id))
    for i=1:length(g.n2o_cell_to_child_id)
      if Gridap.Adaptivity.RefinementRuleType(rrules[i])==Gridap.Adaptivity.WithoutRefinement
        cell_maps_ptrs[i]=wo_refinement_map_id
      else
        cell_maps_ptrs[i]=g.n2o_cell_to_child_id[i]
      end
    end
    return lazy_map(Reindex(cell_maps),cell_maps_ptrs)
  end 
end

const NewToOldReindexArray{T,N}=Gridap.Arrays.LazyArray{<:Fill{<:Reindex},
                                                        T,
                                                        N,
                                                        <:Tuple{<:AbstractVector{<:Integer}}} where {T,N}

function Gridap.Arrays.lazy_map(::typeof(evaluate),
                                f::NewToOldReindexArray{T,N},
                                x::Fill{<:AbstractVector{<:Point}}) where {T,N}
    reindex_map=f.maps.value
    reindex_map_fields=reindex_map.values
    reindex_map_fields_x=[evaluate(reindex_map_fields[i],x.value) for i=1:length(reindex_map_fields)]
    lazy_map(Reindex(reindex_map_fields_x),f.args[1])
end 

function Gridap.Arrays.lazy_map(::typeof(evaluate),
                                f::NewToOldReindexArray{T,N},
                                x::CompressedArray{<:AbstractVector{<:Point}}) where {T,N}
    @notimplementedif length(x.values)>1
    reindex_map=f.maps.value
    reindex_map_fields=reindex_map.values
    reindex_map_fields_x=[evaluate(reindex_map_fields[i],x.values[1]) for i=1:length(reindex_map_fields)]
    lazy_map(Reindex(reindex_map_fields_x),f.args[1])
end

function Gridap.Arrays.lazy_map(::typeof(evaluate),
                                f::Union{Fill{<:AbstractVector{<:Gridap.Fields.Field}},Fill{<:Gridap.Fields.Field}},
                                gx::NewToOldReindexArray{T,N})where {T,N}
    gxmv=gx.maps.value.values
    fv=f.value
    fvgx=[evaluate(fv,gxmv[i]) for i=1:length(gxmv)]
    lazy_map(Reindex(fvgx),gx.args[1])
end

function Gridap.Arrays.lazy_map(::typeof(evaluate),
                                f::Union{CompressedArray{<:AbstractVector{<:Gridap.Fields.Field}},
                                         CompressedArray{<:Gridap.Fields.Field}},
                                gx::NewToOldReindexArray{T,N})where {T,N}
    @notimplementedif length(f.values)>1
    gxmv=gx.maps.value.values
    fv=f.values[1]
    fvgx=[evaluate(fv,gxmv[i]) for i=1:length(gxmv)]
    lazy_map(Reindex(fvgx),gx.args[1])
end

function Gridap.Arrays.lazy_map(m::Gridap.Fields.TransposeMap,
                                gx::NewToOldReindexArray{T,N}) where {T,N}
    gxmv=gx.maps.value.values
    gxmvT=[evaluate(m,gxmv[i]) for i=1:length(gxmv)]
    lazy_map(Reindex(gxmvT),gx.args[1]) 
end


function Gridap.Arrays.lazy_map(m::Gridap.Fields.BroadcastingFieldOpMap,
                                gx::Fill,
                                fx::NewToOldReindexArray{T,N}) where {T,N}
  fxs_op_gx=[evaluate(m,fx.maps.value.values[i],gx.value) for i=1:length(fx.maps.value.values)]
  lazy_map(Reindex(fxs_op_gx),fx.args[1])    
end

function Gridap.Arrays.lazy_map(m::Gridap.Fields.BroadcastingFieldOpMap,
                                gx::CompressedArray,
                                fx::NewToOldReindexArray{T,N}) where {T,N}
  @notimplementedif length(gx.values)>1
  fxs_op_gx=[evaluate(m,fx.maps.value.values[i],gx.values[1]) for i=1:length(fx.maps.value.values)]
  lazy_map(Reindex(fxs_op_gx),fx.args[1])    
end

function Gridap.Arrays.lazy_map(m::Gridap.Fields.BroadcastingFieldOpMap,
                                fx::NewToOldReindexArray{T,N},
                                gx::Fill) where {T,N}
  fxs_op_gx=[evaluate(m,fx.maps.value.values[i],gx.value) for i=1:length(fx.maps.value.values)]
  lazy_map(Reindex(fxs_op_gx),fx.args[1]) 
end

function Gridap.Arrays.lazy_map(m::Gridap.Fields.BroadcastingFieldOpMap,
                                fx::NewToOldReindexArray{T,N},
                                gx::CompressedArray) where {T,N}
  @notimplementedif length(gx.values)>1
  fxs_op_gx=[evaluate(m,fx.maps.value.values[i],gx.values[1]) for i=1:length(fx.maps.value.values)]
  lazy_map(Reindex(fxs_op_gx),fx.args[1])
end

function Gridap.Arrays.lazy_map(m  :: Gridap.Fields.IntegrationMap, 
                                fx :: NewToOldReindexArray{T,N}, 
                                w  :: Fill, 
                                j  :: Fill) where {T,N}
  fxs_op_w_op_j=[evaluate(m,fx.maps.value.values[i],w.value,j.value) for i=1:length(fx.maps.value.values)]
  lazy_map(Reindex(fxs_op_w_op_j),fx.args[1])
end
#### END OPTIMIZATIONS
