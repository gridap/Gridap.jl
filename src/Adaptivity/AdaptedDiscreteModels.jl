
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
function is_child(m1::AdaptedDiscreteModel,m2::DiscreteModel)
  return get_parent(m1) === m2 # m1 = refine(m2)
end

function is_child(m1::AdaptedDiscreteModel,m2::AdaptedDiscreteModel)
  return get_parent(m1) === get_model(m2) # m1 = refine(m2)
end

is_child(m1::DiscreteModel,m2::AdaptedDiscreteModel) = false

is_related(m1::AdaptedDiscreteModel,m2::DiscreteModel) = is_child(m1,m2)
is_related(m1::DiscreteModel,m2::AdaptedDiscreteModel) = is_child(m2,m1)
is_related(m1::AdaptedDiscreteModel,m2::AdaptedDiscreteModel) = is_child(m1,m2) || is_child(m2,m1)

# Model Refining
function refine(model::DiscreteModel,args...;kwargs...) :: AdaptedDiscreteModel
  @abstractmethod
end

function refine(model::AdaptedDiscreteModel,args...;kwargs...)
  ref_model = refine(model.model,args...;kwargs...)
  return AdaptedDiscreteModel(ref_model.model,model,ref_model.glue)
end


# UnstructuredDiscreteModel Refining

function refine(model::UnstructuredDiscreteModel{Dc,Dp}) where {Dc,Dp}
  
  # Create new model
  topo = _refine_unstructured_topology(model.grid_topology)
  grid = UnstructuredGrid(get_vertex_coordinates(topo),get_faces(topo,Dc,0),get_reffes(model.grid),get_cell_type(topo),OrientationStyle(topo))
  nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
  labels = FaceLabeling(nfaces)
  ref_model = UnstructuredDiscreteModel(grid,topo,labels)

  # Create ref glue
  glue = _get_refinement_glue(topo,model.grid_topology)

  return AdaptedDiscreteModel(ref_model,model,glue)
end

_refine_unstructured_topology(topo::UnstructuredGridTopology) = @notimplemented

function _refine_unstructured_topology(topo::UnstructuredGridTopology{Dc,2}) where Dc
  # In dimension D=2, we allow mix and match of TRI and QUAD cells
  @notimplementedif !all(map(p -> p ∈ [QUAD,TRI],topo.polytopes))

  coords_new  = _get_refined_vertex_coordinates(topo)
  c2n_map_new = _get_refined_cell_to_vertex_map(topo)

  nC_old = num_faces(topo,2)
  nC_new = length(c2n_map_new)
  cell_type_new = Vector{Int}(undef,nC_new)
  k = 1
  for iC = 1:nC_old
    p = topo.polytopes[topo.cell_type[iC]]
    range = k:k+num_children(p)-1
    cell_type_new[range] .= topo.cell_type[iC]
    k += num_children(p)
  end
  
  return UnstructuredGridTopology(coords_new,c2n_map_new,cell_type_new,topo.polytopes,topo.orientation_style)
end

function _get_refined_vertex_coordinates(topo::UnstructuredGridTopology{Dc,2}) where Dc
  @notimplementedif !all(map(p -> p ∈ [QUAD,TRI],topo.polytopes))
  coords    = topo.vertex_coordinates
  e2n_map   = get_faces(topo,1,0)
  c2n_map   = get_faces(topo,2,0)

  # Old sizes
  nN_old = num_faces(topo,0) # Nodes
  nE_old = num_faces(topo,1) # Edges

  # Number of Adapted nodes: nN_new = nN + nE + nQUADS
  iQUAD     = findfirst(p -> p == QUAD, topo.polytopes)
  quad_mask = map(i -> i == iQUAD, topo.cell_type)
  nN_new    = nN_old + nE_old + count(quad_mask)

  # Create new coordinates
  coords_new = Vector{eltype(coords)}(undef,nN_new)
  coords_new[1:nN_old] .= coords
  coords_new[nN_old+1:nN_old+nE_old] .= map(pts -> sum(coords[pts])/length(pts), e2n_map)
  coords_new[nN_old+nE_old+1:nN_new] .= map(pts -> sum(coords[pts])/length(pts), c2n_map[findall(quad_mask)])
  return coords_new
end

function _get_refined_cell_to_vertex_map(topo::UnstructuredGridTopology{Dc,2}) where Dc
  @notimplementedif !all(map(p -> p ∈ [QUAD,TRI],topo.polytopes))
  c2n_map = get_faces(topo,2,0)
  c2e_map = get_faces(topo,2,1)

  # Allocate map ptr and data arrays
  nC_new    = sum(ct -> num_children(topo.polytopes[ct]), topo.cell_type)
  nData_new = sum(ct -> num_refined_faces(topo.polytopes[ct]), topo.cell_type)

  ptrs_new  = Vector{Int}(undef,nC_new+1)
  data_new  = Vector{Int}(undef,nData_new)

  # Old sizes
  nN_old = num_faces(topo,0) # Nodes
  nE_old = num_faces(topo,1) # Edges
  nC_old = num_faces(topo,2) # Cells

  k = 1
  ptrs_new[1] = 1
  for iC = 1:nC_old
    p = topo.polytopes[topo.cell_type[iC]]

    # New Node ids from old N,E,C ids
    N = c2n_map[iC]
    E = c2e_map[iC] .+ nN_old
    C = iC + nN_old + nE_old
    new_nodes_ids, p_children = get_refined_children(p,[N,E,C])
    
    nChild = length(new_nodes_ids)
    for iChild = 1:nChild
      ptrs_new[k+1] = ptrs_new[k] + length(new_nodes_ids[iChild])
      data_new[ptrs_new[k]:ptrs_new[k+1]-1] .= new_nodes_ids[iChild]
      k = k+1
    end
  end

  return Table(data_new,ptrs_new)
end

function _get_refinement_glue(ftopo::T,ctopo::T) where {Dc,T<:UnstructuredGridTopology{Dc,2}}
  nC_old = num_faces(ctopo,2)
  nC_new = num_faces(ftopo,2)
  polys  = ctopo.polytopes

  f2c_cell_map      = Vector{Int}(undef,nC_new)
  fcell_to_child_id = Vector{Int}(undef,nC_new)

  k = 1
  for iC = 1:nC_old
    p = polys[ctopo.cell_type[iC]]
    range = k:k+num_children(p)-1
    f2c_cell_map[range] .= iC
    fcell_to_child_id[range] .= collect(1:num_children(p))
    k += num_children(p)
  end

  nP        = length(polys)
  reffes    = map(LagrangianRefFE,Fill(Float64,nP),polys,Fill(1,nP))
  rrules    = map(RefinementRule,reffes,Fill(2,nP))
  fine_to_rrules = lazy_map(Reindex(rrules),ftopo.cell_type)

  f2c_faces_map = [Int[],Int[],f2c_cell_map]
  return AdaptivityGlue(f2c_faces_map,fcell_to_child_id,fine_to_rrules)
end

function num_children(p::Polytope{2})
  (p == QUAD || p == TRI) && (return 4)
  @notimplemented
end

function num_children(p::Polytope{3})
  (p == HEX || p == TET) && (return 8)
  @notimplemented
end

function num_refined_faces(p::Polytope)
  (p == QUAD) && (return num_children(p)*num_faces(p,0))
  (p == TRI)  && (return num_children(p)*num_faces(p,0))
  @notimplemented
end


function get_refined_children(p::Polytope{2}, face_ids)
  N,E,C = face_ids # Node, Edge and Cell global ids

  if p == QUAD 
    ids = [[N[1],E[1],E[3],C],
           [E[1],N[2],C,E[4]],
           [E[3],C,N[3],E[2]],
           [C,E[4],E[2],N[4]]]
    return (ids,QUAD)
  elseif p == TRI
    ids = [[N[1],E[1],E[2]],
           [E[1],N[2],E[3]],
           [E[1],E[2],E[3]],
           [E[2],E[3],N[3]]]
    return (ids,TRI)
  end
  @notimplemented
end

function get_refined_children(p::Polytope{3}, face_ids)
  @notimplemented
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
  rrules         = Fill(RefinementRule(reffe,num_refinements),length(fcell_child_id))
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

