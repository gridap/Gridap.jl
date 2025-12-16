
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

function get_d_to_fface_to_cface(model::AdaptedDiscreteModel)
  ftopo = get_grid_topology(get_model(model))
  ctopo = get_grid_topology(get_parent(model))
  glue  = get_adaptivity_glue(model)
  return get_d_to_fface_to_cface(glue,ctopo,ftopo)
end