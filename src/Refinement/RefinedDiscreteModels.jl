

"""
Refinement glue between two nested triangulations
- `f2c_faces_map`     : Given a fine face gid, returns the gid of the coarse face containing it. 
- `fcell_to_child_id` : Given a fine cell gid, returns the local child id within the coarse cell containing it.
- `f2c_ref_cell_map`  : Ref coordinate map between the fcells and their parent ccell, i.e 
                          for each fine cell defines Φ st. x_c = Φ(x_f)
"""
struct RefinementGlue{A,B,C,D} <: GridapType
  f2c_faces_map           :: A
  fcell_to_child_id       :: B
  f2c_reference_cell_map  :: C
  c2f_faces_map           :: D

  function RefinementGlue(f2c_faces_map,
                          fcell_to_child_id,
                          f2c_reference_cell_map)
    c2f_faces_map = get_c2f_faces_map(f2c_faces_map[end],f2c_reference_cell_map)

    A = typeof(f2c_faces_map)
    B = typeof(fcell_to_child_id)
    C = typeof(f2c_reference_cell_map)
    D = typeof(c2f_faces_map)
    new{A,B,C,D}(f2c_faces_map,fcell_to_child_id,f2c_reference_cell_map,c2f_faces_map)
  end
end

"""
Ref coordinate map from fcell ref coords to ccell ref coords.
Size -> Number of children per coarse cell.
"""
function get_f2c_reference_cell_map(reffe,ref)
  ref_grid = Geometry.UnstructuredGrid(Visualization.compute_reference_grid(reffe,ref))
  return get_cell_map(ref_grid)
end

"""
Ref coordinate map between fine and coarse cells. 
Size -> Number of fine cells. 
"""
function get_f2c_reference_coordinate_map(g::RefinementGlue)
  m = Reindex(g.f2c_reference_cell_map)
  return lazy_map(m,g.fcell_to_child_id)
end

function get_c2f_faces_map(fcell_to_ccell,f2c_reference_cell_map)
  nC = maximum(fcell_to_ccell)
  nF = length(fcell_to_ccell)
  nChildren = length(f2c_reference_cell_map) # TODO: This can be modified for different number of children per cell

  ccell_to_fcell = [fill(-1,nChildren) for i in 1:nC]
  cidx = fill(1,nC)
  for iF in 1:nF
    iC = fcell_to_ccell[iF]
    ccell_to_fcell[iC][cidx[iC]] = iF
    cidx[iC] += 1
  end

  return ccell_to_fcell
end

#################################################################

"""
  RefinedDiscreteModel

  `DiscreteModel` created by refining another `DiscreteModel`. 
  The refinement hierarchy can be traced backwards by following the 
  `parent` pointer chain. This allows the transfer of dofs 
  between `FESpaces` defined on this model and its ancestors using a 
  `ProjectionTransferOperator`.

"""
struct RefinedDiscreteModel{Dc,Dp,A<:DiscreteModel{Dc,Dp},B<:DiscreteModel{Dc,Dp},C<:RefinementGlue} <: DiscreteModel{Dc,Dp}
  model  ::A
  parent ::B
  glue   ::C

  function RefinedDiscreteModel(model::DiscreteModel{Dc,Dp},parent,glue) where {Dc,Dp}
    @check !isa(model,RefinedDiscreteModel)
    A = typeof(model)
    B = typeof(parent)
    C = typeof(glue)
    return new{Dc,Dp,A,B,C}(model,parent,glue)
  end
end

# DiscreteModel API
Geometry.get_grid(model::RefinedDiscreteModel)          = get_grid(model.model)
Geometry.get_grid_topology(model::RefinedDiscreteModel) = get_grid_topology(model.model)
Geometry.get_face_labeling(model::RefinedDiscreteModel) = get_face_labeling(model.model)

# Other getters
get_model(model::RefinedDiscreteModel)  = model.model
get_parent(model::RefinedDiscreteModel{Dc,Dp,A,<:RefinedDiscreteModel}) where {Dc,Dp,A} = get_model(model.parent)
get_parent(model::RefinedDiscreteModel{Dc,Dp,A,B}) where {Dc,Dp,A,B} = model.parent
get_refinement_glue(model::RefinedDiscreteModel) = model.glue

function refine(model::DiscreteModel,args...;kwargs...) :: RefinedDiscreteModel
  @abstractmethod
end

function refine(model::RefinedDiscreteModel,args...;kwargs...)
  ref_model = refine(model.model,args...;kwargs...)
  return RefinedDiscreteModel(ref_model.model,model,ref_model.glue)
end


# UnstructuredDiscreteModelRefining

function refine(model::UnstructuredDiscreteModel{Dc,Dp}) where {Dc,Dp}
  
  # Create new model
  topo = _refine_unstructured_topology(model.grid_topology)
  grid = UnstructuredGrid(get_vertex_coordinates(topo),get_faces(topo,Dc,0),get_reffes(model.grid),get_cell_type(topo),OrientationStyle(topo))
  nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
  labels = FaceLabeling(nfaces)
  ref_model = UnstructuredDiscreteModel(grid,topo,labels)

  # Create ref glue
  glue = _get_refinement_glue(topo,model.grid_topology)

  return RefinedDiscreteModel(ref_model,model,glue)
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

  # Number of refined nodes: nN_new = nN + nE + nQUADS
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

  f2c_cell_map      = Vector{Int}(undef,nC_new)
  fcell_to_child_id = Vector{Int}(undef,nC_new)

  k = 1
  for iC = 1:nC_old
    p = ctopo.polytopes[ctopo.cell_type[iC]]
    range = k:k+num_children(p)-1
    f2c_cell_map[range] .= iC
    fcell_to_child_id[range] .= collect(1:num_children(p))
    k += num_children(p)
  end

  reffe = LagrangianRefFE(QUAD)
  f2c_reference_cell_map = get_f2c_reference_cell_map(reffe,2)

  f2c_faces_map = [Int[],Int[],f2c_cell_map]
  return RefinementGlue(f2c_faces_map,fcell_to_child_id,f2c_reference_cell_map)
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
  ref_cell_map   = get_f2c_reference_cell_map(reffe,ref)
  glue = RefinementGlue(faces_map,fcell_child_id,ref_cell_map)

  # RefinedModel
  return RefinedDiscreteModel(model_ref,model,glue)
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

