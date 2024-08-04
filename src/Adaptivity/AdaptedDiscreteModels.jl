
"""

  `DiscreteModel` created by refining/coarsening another `DiscreteModel`.

  The refinement/coarsening hierarchy can be traced backwards by following the
  `parent` pointer chain. This allows the transfer of dofs
  between `FESpaces` defined on this model and its ancestors.

"""
struct AdaptedDiscreteModel{Dc,Dp,A<:DiscreteModel{Dc,Dp},B<:DiscreteModel{Dc,Dp},C<:AdaptivityGlue} <: DiscreteModel{Dc,Dp}
  model  ::A
  parent ::B
  glue   ::C

  function AdaptedDiscreteModel(model::DiscreteModel{Dc,Dp},parent,glue) where {Dc,Dp}
    @check !isa(model,AdaptedDiscreteModel)
    A = typeof(model)
    B = typeof(parent)
    C = typeof(glue)
    return new{Dc,Dp,A,B,C}(model,parent,glue)
  end
end

# DiscreteModel API
Geometry.get_grid(model::AdaptedDiscreteModel)          = get_grid(model.model)
Geometry.get_grid_topology(model::AdaptedDiscreteModel) = get_grid_topology(model.model)
Geometry.get_face_labeling(model::AdaptedDiscreteModel) = get_face_labeling(model.model)

# Other getters
get_model(model::AdaptedDiscreteModel)  = model.model
get_parent(model::AdaptedDiscreteModel{Dc,Dp,A,<:AdaptedDiscreteModel}) where {Dc,Dp,A} = get_model(model.parent)
get_parent(model::AdaptedDiscreteModel{Dc,Dp,A,B}) where {Dc,Dp,A,B} = model.parent
get_adaptivity_glue(model::AdaptedDiscreteModel) = model.glue

# Relationships
"""
Returns true if m1 is a "child" model of m2, i.e., if m1 is the result of adapting m2
"""
function is_child(m1::AdaptedDiscreteModel,m2::DiscreteModel)
  return get_parent(m1) === m2 # m1 = refine(m2)
end

function is_child(m1::AdaptedDiscreteModel,m2::AdaptedDiscreteModel)
  return get_parent(m1) === get_model(m2) # m1 = refine(m2)
end

is_child(m1::DiscreteModel,m2::AdaptedDiscreteModel) = false

is_related(m1::DiscreteModel,m2::DiscreteModel) = is_child(m1,m2) || is_child(m2,m1)

# Model Adaptation

"""
  function refine(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel

  Returns an `AdaptedDiscreteModel` that is the result of refining the given `DiscreteModel`.
"""
function refine(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel
  @abstractmethod
end

function refine(model::AdaptedDiscreteModel,args...;kwargs...)
  ref_model = refine(model.model,args...;kwargs...)
  return AdaptedDiscreteModel(ref_model.model,model,ref_model.glue)
end

"""
  function coarsen(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel

  Returns an `AdaptedDiscreteModel` that is the result of coarsening the given `DiscreteModel`.
"""
function coarsen(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel
  @abstractmethod
end

"""
  function adapt(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel

  Returns an `AdaptedDiscreteModel` that is the result of adapting (mixed coarsening and refining)
  the given `DiscreteModel`.
"""
function adapt(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel
  @abstractmethod
end

function adapt(model::AdaptedDiscreteModel,args...;kwargs...)
  adapted_model = adapt(model.model,args...;kwargs...)
  return AdaptedDiscreteModel(adapted_model.model,model,adapted_model.glue)
end

# UnstructuredDiscreteModel refining

abstract type AdaptivityMethod end

function refine(model::UnstructuredDiscreteModel,::AdaptivityMethod,args...;kwargs...)
  @abstractmethod
end

# Handle the user's requested choice for refinement
function string_to_refinement(refinement_method::String, model)
  refinement_method == "red_green" && return RedGreenRefinement()
  refinement_method == "nvb" && return NVBRefinement(model)
  refinement_method == "barycentric" && return BarycentricRefinement()
  refinement_method == "simplexify" && return SimplexifyRefinement()
  error("refinement_method $refinement_method not recognized")
end

function refine(model::UnstructuredDiscreteModel,args...;refinement_method="red_green",kwargs...)
  return refine(string_to_refinement(refinement_method, model),model,args...;kwargs...)
end

# CartesianDiscreteModel refining

function refine(model::CartesianDiscreteModel{Dc}, cell_partition::Int=2) where Dc
  partition = Tuple(fill(cell_partition,Dc))
  return refine(model,partition)
end

function refine(model::CartesianDiscreteModel{Dc}, cell_partition::Tuple) where Dc
  desc = Geometry.get_cartesian_descriptor(model)
  nC   = desc.partition

  # Refinement Glue
  f2c_cell_map, fcell_to_child_id = _create_cartesian_f2c_maps(nC,cell_partition)
  faces_map = [(d==Dc) ? f2c_cell_map : Int[] for d in 0:Dc]
  reffe     = LagrangianRefFE(Float64,first(get_polytopes(model)),1)
  rrules    = RefinementRule(reffe,cell_partition)
  glue = AdaptivityGlue(faces_map,fcell_to_child_id,rrules)

  # Refined model
  domain     = _get_cartesian_domain(desc)
  _model_ref = CartesianDiscreteModel(domain,cell_partition.*nC)

  # Propagate face labels
  coarse_labels = get_face_labeling(model)
  coarse_topo   = get_grid_topology(model)
  fine_topo     = get_grid_topology(_model_ref)
  fine_labels   = refine_face_labeling(coarse_labels,glue,coarse_topo,fine_topo)

  model_ref = CartesianDiscreteModel(get_grid(_model_ref),fine_topo,fine_labels)
  return AdaptedDiscreteModel(model_ref,model,glue)
end

function _get_cartesian_domain(desc::CartesianDescriptor{D}) where D
  origin = desc.origin
  corner = origin + VectorValue(desc.sizes .* desc.partition)
  domain = Vector{eltype(origin)}(undef,2*D)
  for d in 1:D
    domain[d*2-1] = origin[d]
    domain[d*2]   = corner[d]
  end
  return Tuple(domain)
end

@generated function _c2v(idx::Union{NTuple{N,T},CartesianIndex{N}},sizes::NTuple{N,T}) where {N,T}
  res = :(idx[1])
  for d in 1:N-1
    ik = :((idx[$(d+1)]-1))
    for k in 1:d
        ik = :($ik * sizes[$k])
    end
    res = :($res + $ik)
  end
  return res
end

@generated function _create_cartesian_f2c_maps(nC::NTuple{N,T},ref::NTuple{N,T}) where {N,T}
  J_f2c   = Meta.parse(prod(["(",["1+(I[$k]-1)÷ref[$k]," for k in 1:N]...,")"]))
  J_child = Meta.parse(prod(["(",["1+(I[$k]-1)%ref[$k]," for k in 1:N]...,")"]))

  return :(begin
    nF = nC .* ref
    f2c_map   = Vector{Int}(undef,prod(nF))
    child_map = Vector{Int}(undef,prod(nF))

    for (i,I) in enumerate(CartesianIndices(nF))
      J_f2c   = $J_f2c
      J_child = $J_child
      f2c_map[i] = _c2v(J_f2c,nC)
      child_map[i] = _c2v(J_child,ref)
    end

    return f2c_map, child_map
  end)
end

function get_d_to_fface_to_cface(model::AdaptedDiscreteModel)
  ftopo = get_grid_topology(get_model(model))
  ctopo = get_grid_topology(get_parent(model))
  glue  = get_adaptivity_glue(model)
  return get_d_to_fface_to_cface(glue,ctopo,ftopo)
end
