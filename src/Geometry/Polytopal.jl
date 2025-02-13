
struct PolytopalGridTopology{Dc,Dp,Tp} <: GridTopology{Dc,Dp}
  vertex_coordinates::Vector{Point{Dp,Tp}}
  n_m_to_nface_to_mfaces::Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}
  polytopes::Vector{GeneralPolytope{Dc,Dp,Tp,Nothing}}
end

function PolytopalGridTopology(
  vertex_coordinates::Vector{<:Point{Dp,T}},
  cell_vertices::Table,
  polytopes::Vector{<:GeneralPolytope{Dc}}
) where {Dc,Dp,T}
  @check length(cell_vertices) == length(polytopes)

  n = Dc+1
  n_m_to_nface_to_mfaces = Matrix{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,n,n)
  n_m_to_nface_to_mfaces[Dc+1,0+1] = cell_vertices
  vertex_cells = generate_cells_around(cell_vertices,length(vertex_coordinates))
  n_m_to_nface_to_mfaces[0+1,Dc+1] = vertex_cells

  PolytopalGridTopology{Dc,Dp,T}(
    vertex_coordinates, n_m_to_nface_to_mfaces, polytopes
  )
end

function PolytopalGridTopology(topo::UnstructuredGridTopology{Dc}) where Dc
  polys, ctypes = get_polytopes(topo), get_cell_type(topo)
  cell_polys = expand_cell_data(polys,ctypes)
  vertex_coords = get_vertex_coordinates(topo)
  polytopes, cell_nodes = generate_polytopes(
    cell_polys,get_faces(topo,Dc,0),vertex_coords
  )
  PolytopalGridTopology(vertex_coords,cell_nodes,polytopes)
end

function oriented_cell_nodes(
  node_coordinates::AbstractVector{<:Point{2}},cell_nodes::Table{T}
) where T
  ptrs = cell_nodes.ptrs
  data = zeros(T, length(cell_nodes.data))
  for cell in 1:length(cell_nodes)
    nodes = view(cell_nodes,cell)
    perm  = orient_nodes(view(node_coordinates,nodes))
    data[ptrs[cell]:ptrs[cell+1]-1] .= nodes[perm]
  end
  return Table(data,ptrs)
end

function generate_polytopes(
  cell_polys::AbstractArray{<:Polytope{2}},cell_nodes,node_coordinates
)
  o_cell_nodes = oriented_cell_nodes(node_coordinates,cell_nodes)
  cell_coordinates = lazy_map(Broadcasting(Reindex(node_coordinates)),o_cell_nodes)
  polytopes = map(Polygon,cell_coordinates)
  return polytopes, o_cell_nodes
end

function generate_polytopes(
  cell_polys::AbstractArray{<:Polytope{3}},cell_nodes,node_coordinates
)
  polytopes = map(cell_polys,cell_nodes) do p, nodes
    GeneralPolytope{3}(node_coordinates[nodes],get_faces(p,1,0))
  end
  return polytopes, cell_nodes
end

function num_faces(g::PolytopalGridTopology{Dc},d::Integer) where Dc
  if d == 0
    length(g.vertex_coordinates)
  elseif d == Dc
    length(g.polytopes)
  else
    face_to_cells = get_faces(g,d,Dc)
    length(face_to_cells)
  end
end

function get_faces(g::PolytopalGridTopology, dimfrom::Integer, dimto::Integer)
  _setup_faces!(g,dimfrom,dimto)
  g.n_m_to_nface_to_mfaces[dimfrom+1,dimto+1]
end

OrientationStyle(::Type{PolytopalGridTopology}) = NonOriented()
get_vertex_coordinates(g::PolytopalGridTopology) = g.vertex_coordinates
get_cell_type(g::PolytopalGridTopology) = Base.OneTo(length(g.polytopes))
get_polytopes(g::PolytopalGridTopology) = g.polytopes

function GridTopology(::Type{<:Polytope{Dc}},topo::PolytopalGridTopology{Dc}) where Dc
  return topo
end

function GridTopology(::Type{<:Polytope{Df}},topo::PolytopalGridTopology{Dc}) where {Df,Dc}
  node_coordinates = get_vertex_coordinates(topo)
  face_to_nodes = get_faces(topo,Df,0)
  if iszero(Df) || isone(Df)
    polytopes = [ifelse(iszero(Df), VERTEX, SEGMENT)]
    face_to_ftype = fill(Int8(1),num_faces(topo,Df))
    ftopo = UnstructuredGridTopology(node_coordinates, face_to_nodes, face_to_ftype, polytopes)
  else
    face_to_polytopes = get_reffaces(Polytope{Df},topo)
    ftopo = PolytopalGridTopology(node_coordinates, face_to_nodes, face_to_polytopes)
  end
  return ftopo
end

############################################################################################

struct PolytopalGrid{Dc,Dp,Tp,Tn} <: Grid{Dc,Dp}
  node_coordinates::Vector{Point{Dp,Tp}}
  cell_node_ids::Table{Int32,Vector{Int32},Vector{Int32}}
  polytopes::Vector{<:GeneralPolytope{Dc,Dp,Tp,Nothing}}
  facet_normal::Tn

  function PolytopalGrid(
    node_coordinates::Vector{Point{Dp,Tp}},
    cell_node_ids::Table{Ti},
    polytopes::Vector{<:GeneralPolytope{Dc}},
    facet_normal=nothing;
  ) where {Dc,Dp,Tp,Ti}
    Tn = typeof(facet_normal)
    new{Dc,Dp,Tp,Tn}(
      node_coordinates, cell_node_ids, polytopes, facet_normal,
    )
  end
end

function PolytopalGrid(grid::Grid{Dc}) where Dc
  node_coordinates = collect1d(get_node_coordinates(grid))
  cell_polys = expand_cell_data(get_polytopes(grid),get_cell_type(grid))
  polytopes, cell_node_ids = generate_polytopes(
    cell_polys,Table(get_cell_node_ids(grid)),node_coordinates
  )
  PolytopalGrid(node_coordinates,cell_node_ids,polytopes)
end

function PolytopalGrid(topo::PolytopalGridTopology{Dc}) where Dc
  PolytopalGrid(get_vertex_coordinates(topo), get_faces(topo,Dc,0), get_polytopes(topo))
end

OrientationStyle(::Type{<:PolytopalGrid}) = NonOriented()
get_reffes(g::PolytopalGrid) = @notimplemented
get_polytopes(g::PolytopalGrid) = g.polytopes
get_cell_type(g::PolytopalGrid) = Base.OneTo(length(g.polytopes))
get_node_coordinates(g::PolytopalGrid) = g.node_coordinates
get_cell_node_ids(g::PolytopalGrid) = g.cell_node_ids
get_cell_map(g::PolytopalGrid) = Fill(GenericField(identity),num_cells(g))

function get_facet_normal(g::PolytopalGrid)
  @assert g.facet_normal != nothing "This Grid does not have information about normals."
  g.facet_normal
end

############################################################################################

struct PolytopalDiscreteModel{Dc,Dp,Tp,Tn} <: DiscreteModel{Dc,Dp}
  grid::PolytopalGrid{Dc,Dp,Tp,Tn}
  grid_topology::PolytopalGridTopology{Dc,Dp,Tp}
  labels::FaceLabeling
end

function PolytopalDiscreteModel(model::DiscreteModel)
  grid = PolytopalGrid(get_grid(model))
  topo = PolytopalGridTopology(get_grid_topology(model))
  labels = get_face_labeling(model)
  return PolytopalDiscreteModel(grid,topo,labels)
end

get_grid(model::PolytopalDiscreteModel) = model.grid
get_grid_topology(model::PolytopalDiscreteModel) = model.grid_topology
get_face_labeling(model::PolytopalDiscreteModel) = model.labels

function compute_face_nodes(model::PolytopalDiscreteModel,d::Integer)
  topo = get_grid_topology(model)
  get_faces(topo,d,0)
end

function compute_face_own_nodes(model::PolytopalDiscreteModel,d::Integer)
  topo = get_grid_topology(model)
  if iszero(d)
    get_faces(topo,d,0)
  else
    empty_table(num_faces(topo,d))
  end
end

function Grid(::Type{ReferenceFE{Df}},model::PolytopalDiscreteModel{Dc}) where {Df,Dc}
  node_coordinates = get_node_coordinates(model)
  face_to_nodes = Table(get_face_nodes(model,Df))
  if iszero(Df) || isone(Df)
    reffes = [ifelse(iszero(Df), VERTEX1, SEG2)]
    face_to_ftype = fill(Int8(1),num_faces(get_grid_topology(model),Df))
    grid = UnstructuredGrid(node_coordinates, face_to_nodes, reffes, face_to_ftype)
  else
    face_to_polytopes = get_reffaces(Polytope{Df},get_grid_topology(model))
    grid = PolytopalGrid(node_coordinates, face_to_nodes, face_to_polytopes)
  end
  return grid
end

function Grid(::Type{ReferenceFE{Dc}},model::PolytopalDiscreteModel{Dc}) where {Dc}
  get_grid(model)
end

############################################################################################
# Voronoi meshes

function orient_nodes(v)
  _angle(c,v) = atan(v[1]-c[1],v[2]-c[2])
  c = mean(v)
  sortperm(v, by = x -> _angle(c,x), rev = true)
end

function voronoi(model::DiscreteModel)
  new_topo = voronoi(get_grid_topology(model))
  new_grid = PolytopalGrid(new_topo)
  new_labels = FaceLabeling(new_topo) # TODO: Map labels from the original model
  return PolytopalDiscreteModel(new_grid,new_topo,new_labels)
end

function voronoi(topo::GridTopology{Dc}) where Dc
  n2c = get_faces(topo,0,Dc)
  n2e = get_faces(topo,0,Dc-1)
  c2n = get_faces(topo,Dc,0)
  e2n = get_faces(topo,Dc-1,0)
  node_coords = get_vertex_coordinates(topo)

  is_bnode = Geometry.get_isboundary_face(topo,0)
  is_bface = Geometry.get_isboundary_face(topo,Dc-1)
  bnodes = findall(is_bnode)
  bfaces = findall(is_bface)
  node_reindex = find_inverse_index_map(bnodes,num_faces(topo,0))    # Not counting the interior nodes for the global dual mesh
  face_reindex = find_inverse_index_map(bfaces,num_faces(topo,Dc-1))  # Not counting the interior faces for the global dual mesh    

  bulk_cen = map(mean,lazy_map(Broadcasting(Reindex(node_coords)), c2n))
  edge_cen = map(mean,lazy_map(Broadcasting(Reindex(node_coords)), e2n))

  new_node_coords = vcat(bulk_cen, edge_cen[bfaces], node_coords[bnodes])
  n_new_cells = num_faces(topo,0)

  ptrs = zeros(Int,n_new_cells+1)
  data = zeros(Int,n_new_cells * 3*2^(Dc-1)) # rough upper-bound

  f_offset = num_faces(topo,Dc)
  n_offset = f_offset + length(bfaces)

  ptrs[1] = 1
  for n in 1:n_new_cells
    if !is_bnode[n]
      # Case 1: Interior node
      #  --> New cell vertices are the centroids of the neighboring cells
      nbr_cells = view(n2c,n)
      ptrs[n+1] = ptrs[n] + length(nbr_cells)
      connectivity = nbr_cells
    else
      # Case 2: Boundary node
      #  --> New cell vertices given by: 
      #      1. The centroids of the neighboring cells
      #      2. The centroids of the neighboring faces
      #      3. The original node
      nbr_cells = view(n2c,n)
      nbr_faces = face_reindex[filter(f -> is_bface[f], view(n2e,n))]
      new_node  = node_reindex[n]
      connectivity = vcat(nbr_cells, nbr_faces .+ f_offset, new_node + n_offset)
    end
    # Orient the vertices counter-clockwise
    new_coords = view(new_node_coords,connectivity)
    perm = orient_nodes(new_coords)
    # Store the connectivity
    ptrs[n+1] = ptrs[n] + length(connectivity)
    data[ptrs[n]:ptrs[n+1]-1] .= connectivity[perm]
  end
  resize!(data,ptrs[end]-1)
  new_c2n = Arrays.Table(data,ptrs)

  # Create polytopes
  Tp = eltype(eltype(new_node_coords))
  polytopes = Vector{GeneralPolytope{Dc,Dc,Tp,Nothing}}(undef,n_new_cells)
  for cell in 1:n_new_cells
    vertices = view(new_c2n,cell)
    polytopes[cell] = Polygon(new_node_coords[vertices])
  end

  return PolytopalGridTopology(new_node_coords,new_c2n,polytopes)
end
