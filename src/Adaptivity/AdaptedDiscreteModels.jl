
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

is_related(m1::DiscreteModel,m2::DiscreteModel) = is_child(m1,m2) || is_parent(m1,m2)

# Model Refining
function refine(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel
  @abstractmethod
end

function refine(model::AdaptedDiscreteModel,args...;kwargs...)
  ref_model = refine(model.model,args...;kwargs...)
  return AdaptedDiscreteModel(ref_model.model,model,ref_model.glue)
end


# UnstructuredDiscreteModel Refining

function refine(model::UnstructuredDiscreteModel{Dc,Dp};cells_to_refine=nothing) where {Dc,Dp}
  # Create new model
  rrules, faces_list = _setup_redgreen_coloring(model.grid_topology,cells_to_refine)
  topo   = _refine_unstructured_topology(model.grid_topology,rrules,faces_list)
  reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  grid   = UnstructuredGrid(get_vertex_coordinates(topo),get_faces(topo,Dc,0),reffes,get_cell_type(topo),OrientationStyle(topo))
  labels = FaceLabeling([num_faces(topo,d) for d in 0:num_cell_dims(topo)])
  ref_model = UnstructuredDiscreteModel(grid,topo,labels)

  # Create ref glue
  glue = _get_refinement_glue(topo,model.grid_topology,rrules)

  return AdaptedDiscreteModel(ref_model,model,glue)
end

_refine_unstructured_topology(topo::UnstructuredGridTopology,rrules,faces_list) = @notimplemented

function _refine_unstructured_topology(topo::UnstructuredGridTopology{Dc,2},
                                       rrules::AbstractVector{<:RefinementRule},
                                       faces_list::Tuple) where Dc
  # In dimension D=2, we allow mix and match of TRI and QUAD cells
  @notimplementedif !all(map(p -> p ∈ [QUAD,TRI],topo.polytopes))

  coords_new  = _get_new_coordinates_from_faces(topo,faces_list)
  c2n_map_new = _get_refined_cell_to_vertex_map(topo,rrules,faces_list)

  nC_old = num_faces(topo,2)
  nC_new = length(c2n_map_new)
  cell_type_new = Vector{Int8}(undef,nC_new)
  polys = unique(lazy_map(rr->first(get_polytopes(rr.ref_grid)), rrules))
  k = 1
  for iC = 1:nC_old
    rr = rrules[iC]
    p  = first(get_polytopes(rr.ref_grid))
    range = k:k+num_subcells(rr)-1
    cell_type_new[range] .= findfirst(x->x==p,polys)
    k += num_subcells(rr)
  end
  
  return UnstructuredGridTopology(coords_new,c2n_map_new,cell_type_new,polys,topo.orientation_style)
end

function _get_refined_cell_to_vertex_map(topo::UnstructuredGridTopology{Dc,2},
                                         rrules::AbstractVector{<:RefinementRule},
                                         faces_list::Tuple) where Dc
  @notimplementedif !all(map(p -> p ∈ [QUAD,TRI],topo.polytopes))
  nN_old,nE_old,nC_old  = num_faces(topo,0),num_faces(topo,1),num_faces(topo,2)
  c2n_map = get_faces(topo,2,0)
  c2e_map = get_faces(topo,2,1)
  ref_edges = faces_list[2]

  # Allocate map ptr and data arrays
  nC_new    = sum(rr -> num_subcells(rr), rrules)
  nData_new = sum(rr -> sum(ids->length(ids),rr.ref_grid.grid.cell_node_ids), rrules)

  ptrs_new  = Vector{Int}(undef,nC_new+1)
  data_new  = Vector{Int}(undef,nData_new)

  if all(lazy_map(rr->RefinementRuleType(rr)<:RedRefinement,rrules))
    edges_reindexing = 1:nE_old
  else
    edges_reindexing = lazy_map(Reindex(find_inverse_index_map(ref_edges,nE_old)),1:nE_old)
  end

  k = 1
  ptrs_new[1] = 1
  C = nN_old + length(ref_edges) + 1
  for iC = 1:nC_old
    rr = rrules[iC]

    # New Node gids from old N,E,C lids
    N = c2n_map[iC]
    E = edges_reindexing[c2e_map[iC]] .+ nN_old
    sub_conn = get_relabeled_connectivity(rr,(N,E,[C]))
    
    nChild = length(sub_conn)
    ptrs_new[k:k+nChild] .=  sub_conn.ptrs .+ (ptrs_new[k] - 1)
    data_new[ptrs_new[k]:ptrs_new[k+nChild]-1] .= sub_conn.data
    k = k+nChild

    _has_interior_point(rr) && (C += 1)
  end

  return Table(data_new,ptrs_new)
end

function _get_refinement_glue(ftopo::T,ctopo::T,rrules) where {Dc,T<:UnstructuredGridTopology{Dc,2}}
  nC_old = num_faces(ctopo,2)
  nC_new = num_faces(ftopo,2)

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

  f2c_faces_map = [Int[],Int[],f2c_cell_map]
  return AdaptivityGlue(f2c_faces_map,fcell_to_child_id,rrules)
end

# Cartesian Mesh refining

function refine(model::CartesianDiscreteModel; num_refinements::Int=2)
  @check num_refinements >= 2
  ref  = num_refinements
  desc = Geometry.get_cartesian_descriptor(model)
  nC   = desc.partition

  @notimplementedif length(nC) != 2
  @notimplementedif any(map(nCi -> nCi != nC[1],nC))

  domain    = _get_cartesian_domain(desc)
  model_ref = CartesianDiscreteModel(domain,ref.*nC)

  # Glue
  faces_map      = [Int[],Int[],_create_f2c_cell_map(nC,ref)]
  fcell_child_id = _create_child_map(nC,ref)
  reffe          = LagrangianRefFE(Float64,QUAD,1)
  rrules         = RefinementRule(reffe,num_refinements)
  glue = AdaptivityGlue(faces_map,fcell_child_id,rrules)

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

function _create_f2c_cell_map(nC::Tuple,ref::Int)
  D  = length(nC)
  nF = nC[1]*ref

  idx = Tuple.(CartesianIndices((nF,nF)))
  a = map((i,j)->(1+(i-1)÷ref,1+(j-1)÷ref),first.(idx),last.(idx))
  b = map((i,j)->(i-1)*nC[1]+j,first.(a),last.(a))
  return Array(reshape(transpose(b),nF*nF))
end

function _create_child_map(nC::Tuple,ref::Int)
  nF = nC[1]*ref
  elem = reshape(collect(1:ref*ref),(ref,ref))
  slice = transpose(repeat(elem,nC[1]))
  mat = repeat(slice,nC[1])
  return Array(reshape(transpose(mat),nF*nF))
end

