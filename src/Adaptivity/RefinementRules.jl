
abstract type RefinementRuleType end
struct GenericRefinement <: RefinementRuleType end
struct WithoutRefinement <: RefinementRuleType end

"""
  Structure representing the map between a single parent cell and its children.

  Contains: 

  - T :: `RefinementRuleType`, indicating the refinement method.
  - poly :: `Polytope`, representing the geometry of the parent cell.
  - ref_grid :: `DiscreteModel` defined on `poly`, giving the parent-to-children cell map. 
"""
struct RefinementRule{P}
  T         :: RefinementRuleType
  poly      :: P
  ref_grid  :: DiscreteModel
  cmaps     :: Vector{<:Field}
  icmaps    :: Vector{<:Field}
  p2c_cache :: Tuple
end

# Note for devs: 
# - The reason why we are saving both the cell maps and the inverse cell maps is to avoid recomputing
#   them when needed. This is needed for performance when the RefinementRule is used for MacroFEs.
#   Also, in the case the ref_grid comes from a CartesianGrid, we save the cell maps as 
#   AffineFields, which are more efficient than the default linear_combinations.
# - We cannot parametrise the RefinementRule by all it's fields, because we will have different types of
#   RefinementRules in a single mesh. It's the same reason why we don't parametrise the ReferenceFE type.

function RefinementRule(
  T::RefinementRuleType,poly::Polytope,ref_grid::DiscreteModel;
  cell_maps=get_cell_map(ref_grid)
)
  ref_trian = Triangulation(ref_grid)
  cmaps = collect1d(cell_maps)
  icmaps = collect1d(lazy_map(Fields.inverse_map,cell_maps))
  # If there are simplices in the mesh, then for a given point, there may exist a symmetric point of it 
  # with respect to the nearest facet, leading to the identification of an incorrect element and resulting
  # in an error. In ther worst case, if this point is the circumcenter, there may be as many wrong points
  # as the number of facets. Therefore we need to find an additional point. 
  polytopes = get_polytopes(ref_grid)
  has_simplex = any(map(is_simplex,polytopes))
  max_num_facets = mapreduce(num_facets,maximum,polytopes)
  num_nearest_vertices = has_simplex ? min(max_num_facets+1,num_nodes(ref_grid)) : 1
  p2c_cache = CellData._point_to_cell_cache(CellData.KDTreeSearch(;num_nearest_vertices),ref_trian)
  return RefinementRule(T,poly,ref_grid,cmaps,icmaps,p2c_cache)
end

function RefinementRule(T::RefinementRuleType,poly::Polytope,ref_grid::Grid)
  ref_model = UnstructuredDiscreteModel(ref_grid)
  cell_maps = get_cell_map(ref_grid)
  return RefinementRule(T,poly,ref_model;cell_maps)
end

function RefinementRule(T::RefinementRuleType,poly::Polytope,ref_grid::CartesianGrid)
  ref_model = UnstructuredDiscreteModel(ref_grid)
  cell_maps = get_cell_map(ref_grid)
  return RefinementRule(T,poly,ref_model;cell_maps)
end

function RefinementRule(reffe::LagrangianRefFE{D},nrefs::Integer;kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(get_polytope(reffe),partition;kwargs...)
end

function RefinementRule(reffe::LagrangianRefFE{D},partition::NTuple{D,Integer};kwargs...) where D
  return RefinementRule(get_polytope(reffe),partition;kwargs...)
end

function RefinementRule(poly::Polytope{D},nrefs::Integer;kwargs...) where D
  partition = tfill(nrefs,Val{D}())
  return RefinementRule(poly,partition;kwargs...)
end

function RefinementRule(poly::Polytope{D},partition::NTuple{D,Integer};kwargs...) where D
  ref_grid = compute_reference_grid(poly,partition)
  return RefinementRule(GenericRefinement(),poly,ref_grid;kwargs...)
end

function Base.show(io::IO,rr::RefinementRule{P}) where P
  T = RefinementRuleType(rr)
  print(io,"RefinementRule{$(rr.poly),$T}")
end

function Base.:(==)(a::RefinementRule,b::RefinementRule)
  A = get_polytope(a) == get_polytope(b)
  B = num_subcells(a) == num_subcells(b)
  C = RefinementRuleType(a) == RefinementRuleType(b)
  return A && B && C
end

ReferenceFEs.get_polytope(rr::RefinementRule) = rr.poly
get_ref_grid(rr::RefinementRule) = rr.ref_grid
num_subcells(rr::RefinementRule) = num_cells(rr.ref_grid)
num_ref_faces(rr::RefinementRule,d::Int) = num_faces(rr.ref_grid,d)
RefinementRuleType(rr::RefinementRule) :: RefinementRuleType = rr.T

Geometry.get_cell_map(rr::RefinementRule) = rr.cmaps
get_inverse_cell_map(rr::RefinementRule) = rr.icmaps

function get_cell_measures(rr::RefinementRule)
  ref_grid  = get_ref_grid(rr)
  ref_trian = Triangulation(UnstructuredDiscreteModel(ref_grid))
  measures  = get_cell_measure(ref_trian)
  M = sum(measures)
  measures /= M
  return measures
end

function Geometry.get_cell_polytopes(rr::Union{RefinementRule,AbstractArray{<:RefinementRule}})
  polys, cell_type = _get_cell_polytopes(rr)
  return CompressedArray(polys,cell_type)
end

function _get_cell_polytopes(rr::RefinementRule)
  ref_grid   = get_ref_grid(rr)
  polys      = get_polytopes(ref_grid)
  cell_types = get_cell_type(ref_grid)
  return polys, cell_types
end

function _get_cell_polytopes(rrules::AbstractArray{<:RefinementRule})
  rr_polys = lazy_map(rr->get_polytopes(rr.ref_grid),rrules)
  rr_cell_type = lazy_map(rr->get_cell_type(rr.ref_grid),rrules)

  # NOTE: The innermost `unique` is to optimize for CompressedArrays
  polys_new = unique(reduce(vcat,unique(rr_polys)))

  # This assumes that new subcells born from a RefinementRule have consecutive gids, such 
  # that the numbering of the new cells is 
  #           gid_new_cell = gid_RefRule_old_cell + child_id_new_cell
  rr2new_cell_type  = lazy_map(vp->map(p->findfirst(x->x==p,polys_new),vp),rr_polys)
  cell_type_new = reduce(vcat,lazy_map((gids,lids)->lazy_map(Reindex(gids),lids),rr2new_cell_type,rr_cell_type))

  return polys_new, cell_type_new
end

# Geometric maps between coarse and fine points

x_to_cell(rr::RefinementRule,x::Point) = CellData._point_to_cell!(rr.p2c_cache,x)

struct CoarseToFinePointMap <: Map end

function Arrays.return_cache(::CoarseToFinePointMap,rr::RefinementRule,x::AbstractVector{<:Point})
  cmaps = get_inverse_cell_map(rr)
  xi_cache = array_cache(x)
  mi_cache = array_cache(cmaps)

  xi = first(x)
  mi = getindex!(mi_cache,cmaps,1)
  zi_cache = Fields.return_cache(mi,xi)
  zi = zero(Fields.return_type(mi,xi))

  T = typeof(zi)
  ptrs = zeros(Int32,num_subcells(rr)+1)
  ids_cache = CachedArray(zeros(Int32,size(x)))
  y_cache = CachedArray(zeros(T,size(x)))
  return xi_cache, mi_cache, zi_cache, y_cache, ids_cache, ptrs
end

function Arrays.evaluate!(cache,::CoarseToFinePointMap,rr::RefinementRule,x::AbstractVector{<:Point})
  xi_cache, mi_cache, zi_cache, y_cache, ids_cache, ptrs = cache
  cmaps = get_inverse_cell_map(rr)
  setsize!(y_cache,size(x))
  setsize!(ids_cache,size(x))
  y, ids = y_cache.array, ids_cache.array

  # First pass: We count the number of points in each subcell, and store the ids
  fill!(ptrs,0)
  for i in eachindex(x)
    xi = getindex!(xi_cache,x,i)
    id = x_to_cell(rr,xi)
    ids[i] = Int32(id)
    ptrs[id+1] += 1
  end
  Arrays.length_to_ptrs!(ptrs)

  # Second pass: We evaluate cmaps on the points in each subcell
  for i in eachindex(x)
    xi = getindex!(xi_cache,x,i)
    id = ids[i]
    mi = getindex!(mi_cache,cmaps,id)
    y[ptrs[id]] = Fields.evaluate!(zi_cache,mi,xi)
    ids[i] = ptrs[id] # Reverse map
    ptrs[id] += 1
  end
  Arrays.rewind_ptrs!(ptrs)

  return Table(y,ptrs), Table(ids,ptrs)
end

# Topological information functions

"
    get_d_to_face_to_child_faces(rr::RefinementRule)

Given a `RefinementRule`, returns for each parent/coarse face the child/fine faces of the 
same dimension that it contains. Therefore, only fine faces at the coarse cell boundary are 
listed in the returned structure.

Returns: [Face dimension][Coarse Face id] -> [Fine faces]
"
function get_d_to_face_to_child_faces(rr::RefinementRule) 
  get_d_to_face_to_child_faces(RefinementRuleType(rr),rr)
end

# Generic version of the function. Spetialisations may exist for some other ref rule types.
# This generic method relies on get_d_to_face_to_parent_face, and simply inverts the map. 
function get_d_to_face_to_child_faces(::RefinementRuleType,rr::RefinementRule)
  d_to_face_to_parent_face, d_to_face_to_parent_face_dim = get_d_to_face_to_parent_face(rr)

  poly = get_polytope(rr)
  Dc   = num_cell_dims(poly)
  d_to_face_to_child_faces = Vector{Vector{Vector{Int32}}}(undef,Dc+1)

  for cface_dim in 0:Dc
    num_cfaces = num_faces(poly,cface_dim)
    cface_to_child_faces = Vector{Vector{Int32}}(undef,num_cfaces)

    parent_faces     = d_to_face_to_parent_face[cface_dim+1]
    parent_faces_dim = d_to_face_to_parent_face_dim[cface_dim+1]
    parent_pairs     = collect(zip(parent_faces,parent_faces_dim))

    for cface in 1:num_cfaces
      cface_to_child_faces[cface] = findall(p -> (p[1] == cface) && (p[2] == cface_dim), parent_pairs)
    end

    d_to_face_to_child_faces[cface_dim+1] = cface_to_child_faces
  end

  return d_to_face_to_child_faces
end

"""
    get_d_to_face_to_parent_face(rr::RefinementRule)

Given a `RefinementRule`, returns for each fine/child face the parent/coarse face
containing it. The parent face can have higher dimension. 

Returns the tuple (A,B) with 
 - A = [Face dimension][Fine Face id] -> [Parent Face]
 - B = [Face dimension][Fine Face id] -> [Parent Face Dimension]
"""
function get_d_to_face_to_parent_face(rr::RefinementRule)
  get_d_to_face_to_parent_face(rr,RefinementRuleType(rr))
end

# Generic version of the function. Spetializations may exist for some other ref rule types.
function get_d_to_face_to_parent_face(rr::RefinementRule,::RefinementRuleType)
  # WARNING: The functions below are NOT general for any point and any polytope.
  #          They are only valid for the specific case of a refinement rule.
  tol = 1.e-10
  function belongs_to_face(::Val{0},::Val{0},fface_coords,cface_coords)
    return norm(fface_coords[1] - cface_coords[1]) < tol
  end
  function belongs_to_face(::Val{0},::Val{1},fface_coords,cface_coords)
    norm(cross(fface_coords[1] - cface_coords[1], fface_coords[1] - cface_coords[2])) < tol
  end
  function belongs_to_face(::Val{0},::Val{2},fface_coords,cface_coords)
    n = cross(cface_coords[2] - cface_coords[1], cface_coords[3] - cface_coords[1])
    return norm(sum(map(ccoords -> dot(n,ccoords - fface_coords[1]), cface_coords))) < tol
  end
  function belongs_to_face(::Val{fface_dim},::Val{cface_dim},fface_coords,cface_coords) where {fface_dim,cface_dim}
    return all(map(p -> belongs_to_face(Val(0),Val(cface_dim),[p],cface_coords),fface_coords))
  end

  ref_grid = get_ref_grid(rr)
  topo = get_grid_topology(ref_grid)
  poly = get_polytope(rr)

  fnode_coords = get_node_coordinates(ref_grid)
  cnode_coords = get_vertex_coordinates(poly)

  Dc = num_cell_dims(ref_grid)
  d_to_face_to_parent_face = Vector{Vector{Int32}}(undef,Dc+1)
  d_to_face_to_parent_face_dim = Vector{Vector{Int32}}(undef,Dc+1)

  # For each fface dimension
  for fface_dim in 0:Dc 
    fface_nodes = Geometry.get_faces(topo,fface_dim,0)
    fface_node_coords = lazy_map(nodes -> lazy_map(Reindex(fnode_coords),nodes),fface_nodes)

    num_ffaces = Geometry.num_faces(topo,fface_dim)
    fface_to_parent_face = fill(Int32(-1),num_ffaces)
    fface_to_parent_face_dim = fill(Int32(-1),num_ffaces)

    # For each fface find the parent face containing it
    for (fface,fcoords) in enumerate(fface_node_coords)
      found = false
      cface_dim = fface_dim
      # Start with cfaces of the same dimension as the fface, and go up until reaching Dc-1
      while (!found) && (cface_dim < Dc)
        cface_nodes = get_faces(poly,cface_dim,0)
        cface_node_coords = lazy_map(nodes -> lazy_map(Reindex(cnode_coords),nodes),cface_nodes)

        for (cface,ccoords) in enumerate(cface_node_coords)
          if !found && belongs_to_face(Val(fface_dim),Val(cface_dim),fcoords,ccoords)
            found = true
            fface_to_parent_face[fface] = cface
            fface_to_parent_face_dim[fface] = cface_dim
          end
        end
        cface_dim += 1
      end
      if !found # Belongs to the cell itself (dimension Dc)
        fface_to_parent_face[fface] = 1
        fface_to_parent_face_dim[fface] = Dc
      end
    end

    d_to_face_to_parent_face[fface_dim+1] = fface_to_parent_face
    d_to_face_to_parent_face_dim[fface_dim+1] = fface_to_parent_face_dim
  end

  return d_to_face_to_parent_face, d_to_face_to_parent_face_dim
end

function _get_terms(poly::Polytope,orders)
  _nodes, facenodes = ReferenceFEs._compute_nodes(poly,orders)
  terms = ReferenceFEs._coords_to_terms(_nodes,orders)
  return terms
end

function _get_face_orders(p::Polytope{Dc},D::Int,orders::Tuple) where Dc
  @check length(orders) == Dc
  @check 1 <= D < Dc
  @check is_n_cube(p)

  if D == 1 # Edges (2D, 3D)
    tangents = get_edge_tangent(p)
    face_orders = map(tangents) do t
      axis = findfirst(i -> abs(t[i]) > 0.5 ,1:Dc)
      return [orders[axis]]
    end
  elseif D == Dc-1 # Faces (3D)
    normals = get_facet_normal(p)
    face_orders = map(normals) do n
      mask = map(i -> abs(n[i]) < 1.e-3,1:Dc)
      return [orders[mask]...]
    end
  else
    @notimplemented
  end

  return face_orders
end

function _get_local_dof_ranges(p::Polytope{Dc},orders) where Dc
  @check length(orders) == Dc
  @check is_n_cube(p)
  idx = CartesianIndices(Tuple(fill(2,Dc)))
  ranges = map(idx) do ii
    map(Tuple(ii),orders) do i,o
      (i-1)*o+1:i*o+1
    end
  end
  return ranges
end

"""
Given a `RefinementRule` of dimension Dc and a Dc-Tuple `fine_orders` of approximation orders, 
returns a map between the fine nodal dofs of order `fine_orders` in the reference grid and the 
coarse nodal dofs of order `2⋅fine_orders` in the coarse parent cell. 

The result is given for each coarse/parent face of dimension `D` as a list of the corresponding 
fine dof lids, i.e 
- [coarse face][coarse dof lid] -> fine dof lid
"""
function get_face_subface_ldof_to_cell_ldof(
  rr::RefinementRule{<:ExtrusionPolytope{Dc}},
  fine_orders::NTuple{Dc,<:Integer},
  D::Int
) where Dc
  poly  = get_polytope(rr)
  coarse_orders = 2 .* fine_orders
  coarse_reffe  = ReferenceFE(poly,lagrangian,Float64,coarse_orders)
  coarse_face_polys = CompressedArray(ReferenceFEs._compute_reffaces_and_face_types(poly,Val(D))...)
  c_edge_to_coarse_dof = coarse_reffe.face_nodes[get_dimranges(poly)[D+1]]
  
  model = get_ref_grid(rr)
  fine_face_grid = Grid(ReferenceFE{D},model)
  fine_face_polys = CompressedArray(get_polytopes(fine_face_grid),get_cell_type(fine_face_grid))

  d_to_face_to_child_faces = get_d_to_face_to_child_faces(rr)
  face_to_child_faces = d_to_face_to_child_faces[D+1]
  
  coarse_face_orders = _get_face_orders(poly,D,coarse_orders)
  fine_face_orders   = _get_face_orders(poly,D,fine_orders)
  
  num_coarse_faces = num_faces(coarse_reffe,D)
  coarse_dofs_above_fine_dofs = Vector{Vector{Vector{Int32}}}(undef,num_coarse_faces)
  for cF in 1:num_coarse_faces
    coarse_face_poly = coarse_face_polys[cF]
    coarse_terms = _get_terms(coarse_face_poly,coarse_face_orders[cF])
    coarse_dofs  = zeros(Int32,Tuple(maximum(coarse_terms)))
    coarse_dofs[coarse_terms] .= c_edge_to_coarse_dof[cF]
    fine_face_to_dof_range = _get_local_dof_ranges(coarse_face_poly,fine_face_orders[cF])
  
    child_faces = face_to_child_faces[cF]
    fine_dofs = Vector{Vector{Int32}}(undef,length(child_faces))
    for (i,fF) in enumerate(child_faces)
      fine_face_poly = fine_face_polys[fF]
      fine_terms = _get_terms(fine_face_poly,fine_face_orders[cF])

      local_dof_range = fine_face_to_dof_range[i]
      local_coarse_dofs = view(coarse_dofs,local_dof_range...)
      fine_dofs[i] = map(Reindex(local_coarse_dofs),fine_terms)
    end
    coarse_dofs_above_fine_dofs[cF] = fine_dofs
  end

  return coarse_dofs_above_fine_dofs
end

"""
    get_cface_to_num_own_ffaces(rr::RefinementRule)

Given a `RefinementRule`, returns for each parent/coarse face the number of child/fine faces of all 
dimensions that it owns. 
"""
function get_cface_to_num_own_ffaces(rr::RefinementRule)
  d_to_face_to_parent_face, d_to_face_to_parent_face_dim = get_d_to_face_to_parent_face(rr)
  return _compute_cface_to_num_own_ffaces(
    rr,d_to_face_to_parent_face,d_to_face_to_parent_face_dim
  )
end

function _compute_cface_to_num_own_ffaces(
  rr::RefinementRule,
  d_to_face_to_parent_face,
  d_to_face_to_parent_face_dim
)
  poly = get_polytope(rr)
  Dc   = num_cell_dims(poly)
  coffsets = get_offsets(poly)

  cface_to_num_own_ffaces = zeros(Int32,num_faces(poly))
  for fface_dim in 0:Dc
    parent_faces     = d_to_face_to_parent_face[fface_dim+1]
    parent_faces_dim = d_to_face_to_parent_face_dim[fface_dim+1]
    for (cface,cface_dim) in zip(parent_faces,parent_faces_dim)
      co = coffsets[cface_dim+1]
      cface_to_num_own_ffaces[cface+co] += 1
    end
  end

  return cface_to_num_own_ffaces
end

"""
    get_cface_to_own_ffaces(rr::RefinementRule)

Given a `RefinementRule`, returns for each parent/coarse face the child/fine faces of all 
dimensions that it owns. 
"""
function get_cface_to_own_ffaces(rr::RefinementRule)
  d_to_face_to_parent_face, d_to_face_to_parent_face_dim = get_d_to_face_to_parent_face(rr)
  return _compute_cface_to_own_ffaces(
    rr,d_to_face_to_parent_face,d_to_face_to_parent_face_dim
  )
end

function _compute_cface_to_own_ffaces(
  rr::RefinementRule,
  d_to_fface_to_cface,
  d_to_fface_to_cface_dim,
  cface_to_num_own_ffaces = _compute_cface_to_num_own_ffaces(
    rr, d_to_fface_to_cface, d_to_fface_to_cface_dim
  )
)
  poly = get_polytope(rr)
  topo = get_grid_topology(rr.ref_grid)
  Dc   = num_cell_dims(poly)

  coffsets = get_offsets(poly)
  foffsets = get_offsets(topo)

  cface_to_own_ffaces = map(nfaces -> zeros(Int32,nfaces),cface_to_num_own_ffaces)
  ptrs = fill(1,num_faces(poly))
  for fface_dim in 0:Dc
    cfaces     = d_to_fface_to_cface[fface_dim+1]
    cfaces_dim = d_to_fface_to_cface_dim[fface_dim+1]
    for (fface,(cface,cface_dim)) in enumerate(zip(cfaces,cfaces_dim))
      co = coffsets[cface_dim+1]
      fo = foffsets[fface_dim+1]
      cface_to_own_ffaces[cface+co][ptrs[cface+co]] = fface+fo
      ptrs[cface+co] += 1
    end
  end

  return cface_to_own_ffaces
end

"""
    get_cface_to_ffaces(rr::RefinementRule)

Given a `RefinementRule`, returns for each parent/coarse face the child/fine faces of all 
dimensions that are on it (owned and not owned).

The implementation aggregates the results of `get_cface_to_own_ffaces`.
"""
function get_cface_to_ffaces(rr::RefinementRule)
  cface_to_own_ffaces = get_cface_to_own_ffaces(rr)
  return aggregate_cface_to_own_fface_data(rr,cface_to_own_ffaces,cface_to_own_ffaces)
end

"""
    aggregate_cface_to_own_fface_data(
        rr::RefinementRule,
        cface_to_own_fface_to_data :: AbstractVector{<:AbstractVector{T}}
    ) where T
  
Given a `RefinementRule`, and a data structure `cface_to_own_fface_to_data` that contains
data for each child/fine face owned by each parent/coarse face, returns a data structure
`cface_to_fface_to_data` that contains the data for each child/fine face contained in the 
closure of each parent/coarse face (i.e the fine faces are owned and not owned).

The implementation makes sure that the resulting data structure is ordered according to the
fine face numbering in `get_cface_to_ffaces(rrule)` (which in turn is by increasing fine face id).
"""
function aggregate_cface_to_own_fface_data(
  rr::RefinementRule,
  cface_to_own_fface_to_data :: AbstractVector{<:AbstractVector{T}}
) where T
  cface_to_own_ffaces = get_cface_to_own_ffaces(rr)
  return aggregate_cface_to_own_fface_data(rr,cface_to_own_ffaces,cface_to_own_fface_to_data)
end

function aggregate_cface_to_own_fface_data(
  rr::RefinementRule,
  cface_to_own_ffaces::AbstractVector{Vector{Int32}},
  cface_to_own_fface_to_data :: AbstractVector{<:AbstractVector{T}}
) where T
  poly = get_polytope(rr)
  Dc   = num_cell_dims(poly)
  coffsets = get_offsets(poly)

  # For each cface, we concatenate the entries of face_own_dfaces of any 
  # other cface of lower (or equal) dimension that is owned by the cface
  cface_to_fface = [Int32[] for cface in 1:num_faces(poly)]
  cface_to_fface_to_data = [T[] for cface in 1:num_faces(poly)]
  for d1 in 0:Dc
    o1 = coffsets[d1+1]
    for d2 in 0:d1
      o2 = coffsets[d2+1]
      d1_to_d2_faces = get_faces(poly,d1,d2)
      for (d1_cface,d2_cfaces) in enumerate(d1_to_d2_faces)
        append!(
          cface_to_fface[d1_cface + o1],
          cface_to_own_ffaces[d2_cfaces .+ o2]...
        )
        append!(
          cface_to_fface_to_data[d1_cface + o1],
          cface_to_own_fface_to_data[d2_cfaces .+ o2]...
        )
      end
    end
  end

  # We now need to sort the entries of cface_to_fface_to_data according to the order of the
  # fine faces numbering in cface_to_fface.
  perms = map(sortperm,cface_to_fface)
  map(permute!,cface_to_fface_to_data,perms)
  return cface_to_fface_to_data
end

"""
    get_cface_to_ffaces_to_lnodes(rr::RefinementRule)

Given a `RefinementRule`, returns

  [coarse face][local child face] -> local fine node ids

where local refers to the local fine numbering within the coarse face.
"""
function get_cface_to_ffaces_to_lnodes(rrule::RefinementRule)
  topo = get_grid_topology(rrule.ref_grid)
  cface_to_ffaces = get_cface_to_ffaces(rrule)
  cface_to_fnodes = get_face_vertices(rrule)
  fface_to_fnodes = get_face_vertices(topo)

  return _compute_cface_to_ffaces_to_lnodes(
    cface_to_ffaces,cface_to_fnodes,fface_to_fnodes
  )
end

"""
    get_cface_to_own_ffaces_to_lnodes(rr::RefinementRule)

Given a `RefinementRule`, returns

  [coarse face][local owned child face] -> local fine node ids
  
where local refers to the local fine numbering within the coarse face.
"""
function get_cface_to_own_ffaces_to_lnodes(rrule::RefinementRule)
  topo = get_grid_topology(rrule.ref_grid)
  cface_to_ffaces = get_cface_to_own_ffaces(rrule)
  cface_to_fnodes = get_face_vertices(rrule)
  fface_to_fnodes = get_face_vertices(topo)

  return _compute_cface_to_ffaces_to_lnodes(
    cface_to_ffaces,cface_to_fnodes,fface_to_fnodes
  )
end

function _compute_cface_to_ffaces_to_lnodes(
  cface_to_ffaces::AbstractVector{Vector{Int32}},
  cface_to_fnodes::AbstractVector{Vector{Int32}},
  fface_to_fnodes::AbstractVector{Vector{Int32}},
)
  ncfaces = length(cface_to_ffaces)
  cface_to_ffaces_to_lnodes = Vector{Vector{Vector{Int32}}}(undef,ncfaces)
  for cface in 1:ncfaces
    ffaces = cface_to_ffaces[cface]
    cface_fnodes = cface_to_fnodes[cface]
    @assert issorted(cface_fnodes)

    ffaces_to_lnodes = Vector{Vector{Int32}}(undef,length(ffaces))
    for (i,fface) in enumerate(ffaces)
      fface_fnodes = fface_to_fnodes[fface]
      lnodes = zeros(Int32, length(fface_fnodes))
      for (j,fnode) in enumerate(fface_fnodes)
        lnode = searchsortedfirst(cface_fnodes,fnode)
        @assert lnode != length(cface_fnodes)+1
        lnodes[j] = lnode
      end
      ffaces_to_lnodes[i] = lnodes
    end
    cface_to_ffaces_to_lnodes[cface] = ffaces_to_lnodes
  end

  return cface_to_ffaces_to_lnodes
end

"""
    ReferenceFEs.get_face_vertices(rr::RefinementRule)

Given a `RefinementRule`, returns for each parent/coarse face the ids of the
child/fine vertices it contains.
"""
function ReferenceFEs.get_face_vertices(rrule::RefinementRule)
  cface_to_ffaces = get_cface_to_ffaces(rrule)
  topo = get_grid_topology(rrule.ref_grid)
  node_range = get_dimrange(topo,0)
  face_vertices = map(ffaces -> filter(fface -> fface ∈ node_range,ffaces), cface_to_ffaces)
  return face_vertices
end

"""
    ReferenceFEs.get_face_coordinates(rr::RefinementRule)

Given a `RefinementRule`, returns for each parent/coarse face the coordinates of the
child/fine vertices it contains.
"""
function ReferenceFEs.get_face_coordinates(rrule::RefinementRule)
  face_vertices = get_face_vertices(rrule)
  topo = get_grid_topology(rrule.ref_grid)
  coords = get_vertex_coordinates(topo)
  face_coords = lazy_map(Broadcasting(Reindex(coords)), face_vertices)
  return face_coords
end

"""
    ReferenceFEs.get_vertex_permutations(rr::RefinementRule)

Given a `RefinementRule`, returns all possible permutations of the child/fine vertices
within the cell.
"""
function ReferenceFEs.get_vertex_permutations(rrule::RefinementRule)
  poly = get_polytope(rrule)
  topo = get_grid_topology(rrule.ref_grid)

  cnodes  = collect(1:num_faces(poly,0))
  fnodes  = collect(1:num_faces(topo,0))
  fcoords = get_node_coordinates(rrule.ref_grid)
  ccoords = get_vertex_coordinates(poly)
  pindex_to_cnodes = get_vertex_permutations(poly)

  pindex_to_fnodes = _compute_face_vertex_permutations(
    poly,cnodes,fnodes,ccoords,fcoords,pindex_to_cnodes
  )
  return pindex_to_fnodes
end

"""
    ReferenceFEs.get_face_vertex_permutations(rr::RefinementRule)

Given a `RefinementRule`, returns for each parent/coarse face the possible permutations of the
child/fine vertices it contains.
"""
function ReferenceFEs.get_face_vertex_permutations(rrule::RefinementRule)
  poly = get_polytope(rrule)
  cface_to_cnodes = get_face_vertices(poly)
  cface_to_fnodes = get_face_vertices(rrule)
  cface_to_ccoords = get_face_coordinates(poly)
  cface_to_fcoords = get_face_coordinates(rrule)
  cface_to_pindex_to_cnodes = get_face_vertex_permutations(poly)

  Ti = eltype(first(first(cface_to_pindex_to_cnodes)))
  cface_to_pindex_to_fnodes = Vector{Vector{Vector{Ti}}}(undef,num_faces(poly))
  for cface in 1:num_faces(poly)
    cface_to_pindex_to_fnodes[cface] = _compute_face_vertex_permutations(
      poly,
      cface_to_cnodes[cface],
      cface_to_fnodes[cface],
      cface_to_ccoords[cface],
      cface_to_fcoords[cface],
      cface_to_pindex_to_cnodes[cface]
    )
  end
  
  return cface_to_pindex_to_fnodes
end

# Implementation comment (Jordi): 
# The function below computes the permutations of the fine nodes on a face, given a permutation 
# of the coarse nodes on the same face. This is done by comparing coordinates. In general, this is 
# not the best idea, since we are doing float comparisons. However, I believe this should not be a problem
# in this case, since RefinementRules should never have two nodes that are too close to each other.
# I have also not found another way of doing this.
function _compute_face_vertex_permutations(
  poly   :: Polytope,                    # The RefinementRule polytope (NOT the face polytope)
  cnodes :: Vector{<:Integer},           # The coarse node ids on the coarse face
  fnodes :: Vector{<:Integer},           # The fine node ids on the coarse face
  ccoords::Vector{VectorValue{D,T}},     # The coarse node coordinates (for the cnodes)
  fcoords::Vector{VectorValue{D,T}},     # The fine node coordinates (for the fnodes)
  pindex_to_cnodes::Vector{Vector{Ti}},  # The coarse node permutations (in local face numbering)
) where {D,T,Ti}
  # Collect the D-dimensional shape functions for the coarse face
  reffe = LagrangianRefFE(Float64,poly,1)
  shapefuns = get_shapefuns(reffe)[cnodes]

  # For each permutation of the coarse face nodes
  pindex_to_fnodes = Vector{Vector{Ti}}(undef,length(pindex_to_cnodes))
  for (pindex, p_cnodes) in enumerate(pindex_to_cnodes)
    # Compute the coordinates of the permuted fine nodes
    p_ccoords = ccoords[p_cnodes]
    cmap = linear_combination(p_ccoords,shapefuns)
    p_fcoords = evaluate(cmap, fcoords)

    # Match the permuted coordinates to find the fine node permutation
    p_fnodes = zeros(Ti,length(fcoords))
    for (i,c) in enumerate(p_fcoords)
      p = findfirst(x -> norm(x-c) < eps(T),fcoords)
      @assert !isnothing(p)
      p_fnodes[i] = p # for non-local numbering, we would do fnodes[p]
    end
    @check !any(iszero,p_fnodes) "Face vertex permutation not found!"

    pindex_to_fnodes[pindex] = p_fnodes
  end
  
  return pindex_to_fnodes
end

"""
    get_cface_to_fface_permutations(rrule::RefinementRule)
    get_cface_to_own_fface_permutations(rrule::RefinementRule)

Given a `RefinementRule`, this function returns: 

  - `cface_to_cpindex_to_ffaces` : For each coarse face, for each coarse face permutation, 
    the permuted ids of the fine faces.
  - `cface_to_cpindex_to_fpindex` : For each coarse face, for each coarse face permutation, 
    the sub-permutation of the fine faces.

The idea is the following: A permutation on a coarse face induces a 2-level permutation for 
the fine faces, i.e 

  - First, the fine faces get shuffled amongs themselves.
  - Second, each fine face has it's orientation changed (given by a sub-permutation).

For instance, let's consider a 1D example, where a SEGMENT is refined into 2 segments: 

       3             4     5
    X-----X  -->  X-----X-----X 
    1     2       1     3     2

Then when aplying the coarse node permutation (1,2) -> (2,1), 
we get the following fine face permutation: 

  - Faces (1,2,3,4,5) get mapped to (2,1,3,5,4)
  - Moreover, the orientation of faces 3 and 5 is changed, 
    i.e we get the sub-permutation (1,1,1,2,2)
"""
function get_cface_to_fface_permutations(rrule::RefinementRule)
  cface_to_ffaces = get_cface_to_ffaces(rrule)
  cface_to_fface_to_clnodes = get_cface_to_ffaces_to_lnodes(rrule)
  return _compute_cface_to_fface_permutations(
    rrule,cface_to_ffaces,cface_to_fface_to_clnodes
  )
end

function get_cface_to_own_fface_permutations(rrule::RefinementRule)
  cface_to_ffaces = get_cface_to_own_ffaces(rrule)
  cface_to_fface_to_clnodes = get_cface_to_own_ffaces_to_lnodes(rrule)
  return _compute_cface_to_fface_permutations(
    rrule,cface_to_ffaces,cface_to_fface_to_clnodes
  )
end

function _compute_cface_to_fface_permutations(
  rrule::RefinementRule,
  cface_to_ffaces,
  cface_to_fface_to_clnodes
)
  poly = get_polytope(rrule)
  topo = get_grid_topology(rrule.ref_grid)
  cface_to_cpindex_to_clnodes = get_face_vertex_permutations(rrule)

  fface_to_polytope = CompressedArray(get_reffaces(topo),get_face_type(topo))
  fface_to_fpindex_to_flnodes = lazy_map(get_vertex_permutations,fface_to_polytope)

  # For each coarse face
  cface_to_cpindex_to_ffaces = Vector{Vector{Vector{Int32}}}(undef,num_faces(poly))
  cface_to_cpindex_to_fpindex = Vector{Vector{Vector{Int32}}}(undef,num_faces(poly))
  for cface in 1:num_faces(poly)
    cpindex_to_clnodes = cface_to_cpindex_to_clnodes[cface]
    fface_to_clnodes = cface_to_fface_to_clnodes[cface]
    flface_to_fpindex_to_flnodes = view(fface_to_fpindex_to_flnodes,cface_to_ffaces[cface])

    # For each coarse face permutation
    cpindex_to_ffaces = Vector{Vector{Int32}}(undef,length(cpindex_to_clnodes))
    cpindex_to_fpindex = Vector{Vector{Int32}}(undef,length(cpindex_to_clnodes))
    for (cpindex,clnodes) in enumerate(cpindex_to_clnodes)
      permuted_ffaces = Vector{Int32}(undef,length(fface_to_clnodes))
      permuted_fpindex = Vector{Int32}(undef,length(fface_to_clnodes))
      # For each fine face within the coarse face
      for (fface,unpermuted_clnodes) in enumerate(fface_to_clnodes)
        # Find the permuted position of the fface within the cface:
        # To do this, we match the permuted nodes of the selected fine face 
        # with the unpermuted nodes of all fine faces within the coarse face.
        permuted_clnodes = cpindex_to_clnodes[cpindex][unpermuted_clnodes]
        permuted_fface = findfirst(
          clnodes -> all(map(pnode -> pnode ∈ clnodes, permuted_clnodes)),
          fface_to_clnodes
        )
        @assert !isnothing(permuted_fface)
        # Find the inner fface permutation index:
        # To do this, we match the unpermuted nodes of the new fine face position 
        # with the permutations of the original fine face.
        permuted_fface_unpermuted_clones = fface_to_clnodes[permuted_fface]
        fpindex_to_flnodes = flface_to_fpindex_to_flnodes[fface]
        fpindex = findfirst(
          flnodes -> all(permuted_fface_unpermuted_clones[flnodes] .== permuted_clnodes),
          fpindex_to_flnodes
        )
        @assert !isnothing(fpindex)
        permuted_ffaces[fface] = permuted_fface
        permuted_fpindex[fface] = fpindex
      end
      cpindex_to_ffaces[cpindex] = permuted_ffaces
      cpindex_to_fpindex[cpindex] = permuted_fpindex
    end
    cface_to_cpindex_to_ffaces[cface] = cpindex_to_ffaces
    cface_to_cpindex_to_fpindex[cface] = cpindex_to_fpindex
  end

  return cface_to_cpindex_to_ffaces,cface_to_cpindex_to_fpindex
end

# Tests 

function test_refinement_rule(rr::RefinementRule; debug=false)
  poly = get_polytope(rr)
  D = num_dims(poly)

  cmaps = get_cell_map(rr)
  inv_cmaps = Adaptivity.get_inverse_cell_map(rr)

  pts = map(x -> VectorValue(rand(D)),1:10)
  if Geometry.is_simplex(poly)
    filter!(p -> sum(p) < 1.0, pts)
  end

  # Checking that inv_cmaps ∘ cmaps = identity
  for p in pts
    ichild = Adaptivity.x_to_cell(rr,p)
    m = cmaps[ichild]
    m_inv = inv_cmaps[ichild]

    y = evaluate(m,p)
    z = evaluate(m_inv,y)
    debug && println(ichild, " :: ", p," -> ",y, " -> ", z, " - ", p ≈ z)
    @test norm(p-z) < 1.e-6
  end

  cell_measures = get_cell_measures(rr)
  cell_polys = Geometry.get_cell_polytopes(rr)

  return nothing
end