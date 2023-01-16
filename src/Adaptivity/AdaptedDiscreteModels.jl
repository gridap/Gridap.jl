
"""
  AdaptedDiscreteModel

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
function refine(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel
  @abstractmethod
end

function refine(model::AdaptedDiscreteModel,args...;kwargs...)
  ref_model = refine(model.model,args...;kwargs...)
  return AdaptedDiscreteModel(ref_model.model,model,ref_model.glue)
end

function coarsen(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel
  @abstractmethod
end

function adapt(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel
  @abstractmethod
end

function adapt(model::AdaptedDiscreteModel,args...;kwargs...)
  adapted_model = adapt(model.model,args...;kwargs...)
  return AdaptedDiscreteModel(adapted_model.model,model,adapted_model.glue)
end

# UnstructuredDiscreteModel Refining

function refine(model::UnstructuredDiscreteModel{Dc,Dp};cells_to_refine=nothing) where {Dc,Dp}
  # Create new model
  rrules, faces_list = setup_edge_based_rrules(model.grid_topology,cells_to_refine)
  topo   = refine_unstructured_topology(model.grid_topology,rrules,faces_list)
  reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  grid   = UnstructuredGrid(get_vertex_coordinates(topo),get_faces(topo,Dc,0),reffes,get_cell_type(topo),OrientationStyle(topo))
  labels = FaceLabeling(topo)
  ref_model = UnstructuredDiscreteModel(grid,topo,labels)

  # Create ref glue
  glue = get_refinement_glue(topo,model.grid_topology,rrules)

  return AdaptedDiscreteModel(ref_model,model,glue)
end

function refine_unstructured_topology(topo::UnstructuredGridTopology{Dc},
                                      rrules::AbstractVector{<:RefinementRule},
                                      faces_list::Tuple) where {Dc}
  coords_new  = get_new_coordinates_from_faces(topo,faces_list)
  c2n_map_new = get_refined_cell_to_vertex_map(topo,rrules,faces_list)
  polys_new, cell_type_new = get_refined_polytopes(rrules)

  return UnstructuredGridTopology(coords_new,c2n_map_new,cell_type_new,polys_new,topo.orientation_style)
end

function get_refined_polytopes(rrules::AbstractArray{<:RefinementRule})
  rr_polys = lazy_map(rr->get_polytopes(rr.ref_grid),rrules)

  # NOTE: The innermost `unique` is to optimize for CompressedArrays
  polys_new = unique(reduce(vcat,unique(rr_polys)))

  rr_cell_type = lazy_map(rr->get_cell_type(rr.ref_grid),rrules)
  rr2new_cell_type  = lazy_map(vp->map(p->findfirst(x->x==p,polys_new),vp),rr_polys)
  cell_type_new = reduce(vcat,lazy_map((gids,lids)->lazy_map(Reindex(gids),lids),rr2new_cell_type,rr_cell_type))

  return polys_new, cell_type_new
end

function get_refinement_glue(ftopo::T,ctopo::T,rrules::AbstractVector{<:RefinementRule}) where {Dc,T<:UnstructuredGridTopology{Dc}}
  nC_old = num_faces(ctopo,Dc)
  nC_new = num_faces(ftopo,Dc)

  f2c_cell_map      = Vector{Int}(undef,nC_new)
  fcell_to_child_id = Vector{Int}(undef,nC_new)

  k = 1
  for iC = 1:nC_old
    rr = rrules[iC]
    range = k:k+num_subcells(rr)-1
    f2c_cell_map[range] .= iC
    fcell_to_child_id[range] .= collect(1:num_subcells(rr))
    k += num_subcells(rr)
  end

  f2c_faces_map = [(d==Dc) ? f2c_cell_map : Int[] for d in 0:Dc]
  return AdaptivityGlue(f2c_faces_map,fcell_to_child_id,rrules)
end

# Cartesian Mesh refining

function refine(model::CartesianDiscreteModel{Dc}, cell_partition::Int=2) where Dc
  partition = Tuple(fill(cell_partition,Dc))
  return refine(model,partition)
end

function refine(model::CartesianDiscreteModel{Dc}, cell_partition::Tuple) where Dc
  desc = Geometry.get_cartesian_descriptor(model)
  nC   = desc.partition

  # Refined model
  domain    = _get_cartesian_domain(desc)
  model_ref = CartesianDiscreteModel(domain,cell_partition.*nC)

  # Glue
  f2c_cell_map, fcell_to_child_id = _create_cartesian_f2c_maps(nC,cell_partition)
  faces_map      = [(d==Dc) ? f2c_cell_map : Int[] for d in 0:Dc]
  reffe          = LagrangianRefFE(Float64,first(get_polytopes(model)),1)
  rrules         = RefinementRule(reffe,cell_partition)
  glue = AdaptivityGlue(faces_map,fcell_to_child_id,rrules)

  return AdaptedDiscreteModel(model_ref,model,glue)
end

function _get_cartesian_domain(desc::CartesianDescriptor{D}) where D
  origin = desc.origin
  corner = origin + VectorValue(desc.sizes .* desc.partition)
  domain = Vector{Int}(undef,2*D)
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