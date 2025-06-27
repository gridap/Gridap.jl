# Two different plot functions depending on the need of having a continuous or discontinuous mesh.
function to_plot_dg_mesh(grid::Grid)
  UnstructuredGrid(grid) |> to_plot_dg_mesh
end

function to_plot_dg_mesh(grid::UnstructuredGrid)
    to_simplex_grid(grid) |> to_dg_mesh
end

function to_plot_mesh(grid::Grid)
  UnstructuredGrid(grid) |> to_plot_mesh
end

function to_plot_mesh(grid::UnstructuredGrid)
    to_simplex_grid(grid) |> to_mesh
end

# For the moment, we simplexify all grids (missing support for quads or general polytopes).
function to_simplex_grid(grid::UnstructuredGrid)
  reffes = get_reffes(grid)
  polys = map(get_polytope,reffes)
  all(is_simplex,polys) ? grid : simplexify(grid)
end

# Create continuous GeometryBasics mesh equivalent of our grid.
function to_mesh(grid::UnstructuredGrid)
  reffes = get_reffes(grid)
  @assert all(reffe->get_order(reffe)==1,reffes)
  @assert all(reffe->is_simplex(get_polytope(reffe)),reffes)
  xs = get_node_coordinates(grid)
  Tp = eltype(eltype(xs))
  Dp = length(eltype(xs))
  ps = collect(reinterpret(GeometryBasics.Point{Dp,Tp},xs))
  cns = get_cell_node_ids(grid)
  Tc = eltype(eltype(cns))
  Dc = num_cell_dims(grid)
  fs = collect(lazy_map(GeometryBasics.NgonFace{Dc+1,Tc},cns))
  GeometryBasics.Mesh(ps,fs)
end

# Create discontinuous grid with the corresponding handling of nodal or cell fields.
function to_point(x::VectorValue{D,T}) where {D,T}
  GeometryBasics.Point{D,T}(Tuple(x))
end

function to_dg_points(grid::UnstructuredGrid)
  node_x = get_node_coordinates(grid)
  coords = to_dg_node_values(grid, node_x)
  [to_point(x) for x in coords]
end

function to_dg_mesh(grid::UnstructuredGrid)
  ps = to_dg_points(grid)
  Dc = num_cell_dims(grid)
  Tc = Int32
  cns = Vector{Vector{Tc}}(undef, num_cells(grid))
  i = 1
  for cell in 1:num_cells(grid)
    cn = zeros(Tc, Dc+1)
    for lnode in 1:(Dc+1)
      cn[lnode] = i
      i += one(Tc)
    end
    cns[cell] = cn
  end
  fs = lazy_map(GeometryBasics.NgonFace{Dc+1,Tc}, cns) |> collect
  GeometryBasics.Mesh(ps,fs)
end

function to_dg_node_values(grid::Grid, node_value::Vector)
  cell_nodes = get_cell_node_ids(grid)
  i = 0
  cache = array_cache(cell_nodes)
  for cell in 1:length(cell_nodes)
    nodes = getindex!(cache, cell_nodes, cell)
    for node in nodes
      i += 1
    end
  end
  T = eltype(node_value)
  values = zeros(T,i)
  i = 0
  for cell in 1:length(cell_nodes)
    nodes = getindex!(cache, cell_nodes, cell)
    for node in nodes
      i += 1
      values[i] = node_value[node]
    end
  end
  values
end

function to_dg_cell_values(grid::Grid, cell_value::Vector)
  cell_nodes = get_cell_node_ids(grid)
  i = 0
  cache = array_cache(cell_nodes)
  for cell in 1:length(cell_nodes)
    nodes = getindex!(cache, cell_nodes, cell)
    for node in nodes
      i += 1
    end
  end
  T = eltype(cell_value)
  values = zeros(T,i)
  i = 0
  for cell in 1:length(cell_nodes)
    nodes = getindex!(cache, cell_nodes, cell)
    for node in nodes
      i += 1
      values[i] = cell_value[cell]
    end
  end
  values
end

# Obtain edge and vertex skeletons:
function to_lowdim_grid(grid::Grid, ::Val{D}) where D
  topo = GridTopology(grid)
  labels = FaceLabeling(topo)
  model = DiscreteModel(grid, topo, labels)
  Grid(ReferenceFE{D}, model)
end

to_lowdim_grid(grid::Grid{D}, ::Val{D}) where D = grid

to_face_grid(grid::Grid)   = to_lowdim_grid(grid, Val(2))
to_edge_grid(grid::Grid)   = to_lowdim_grid(grid, Val(1))
to_vertex_grid(grid::Grid) = to_lowdim_grid(grid, Val(0))

function to_face_grid_with_map(grid::Grid)
  function _fixnodeids!(face_to_nodes,face_to_m,face_to_n)
    for face in 1:length(face_to_nodes)
      m = face_to_m[face]
      n = face_to_n[face]
      p = face_to_nodes.ptrs[face]
      if m⋅n < 0
        a = face_to_nodes.data[p]
        face_to_nodes.data[p] = face_to_nodes.data[p+1]
        face_to_nodes.data[p+1] = a
      end
    end
  end
  D = 2
  topo = GridTopology(grid)
  labels = FaceLabeling(topo)
  model = DiscreteModel(grid, topo, labels)
  nfaces = num_faces(model,D)
  Γ = BoundaryTriangulation(model,IdentityVector(nfaces))
  n_Γ = get_normal_vector(Γ)
  x_Γ = CellPoint(Fill(VectorValue(0.0,0.0),nfaces),Γ,ReferenceDomain())
  face_to_n = collect(n_Γ(x_Γ))
  face_to_ϕ = get_cell_map(Γ)
  face_to_Jt = lazy_map(∇,face_to_ϕ)
  Jt = GenericCellField(face_to_Jt,Γ,ReferenceDomain())
  m_Γ = Operation(a->VectorValue(a[1],a[3],a[5])×VectorValue(a[2],a[4],a[6]))(Jt)
  face_to_m = collect(m_Γ(x_Γ))
  g = UnstructuredGrid(Grid(ReferenceFE{D}, model))
  face_to_nodes = get_cell_node_ids(g)
  _fixnodeids!(face_to_nodes,face_to_m,face_to_n)
  face_to_cells = get_faces(topo,D,num_cell_dims(grid))
  face_to_cell = lazy_map(getindex,face_to_cells,Fill(1,nfaces))
  g,face_to_cell
end

# Obtain grid and cellfield from triangulation:
function to_grid(trian::Triangulation)
  vds = visualization_data(trian,"")
  first(vds).grid
end

function to_grid(trian::Triangulation, uh::Any)
  vds = visualization_data(trian, "", cellfields=[""=>uh])
  grid = first(vds).grid
  nodaldata = first(first(vds).nodaldata)[2]
  scalarnodaldata = map(to_scalar, nodaldata)
  grid, scalarnodaldata
end

to_scalar(x) = x
to_scalar(x::VectorValue) = norm(x)

function to_point1D(trian::Triangulation{Dc,1}, uh::Any) where Dc
  vds = first(visualization_data(trian, "", cellfields=["uh"=>uh]))
  y = vds.nodaldata["uh"]
  x = getindex.(get_node_coordinates(vds.grid),1)
  x, y
end

function to_point1D(trian::Triangulation{Dc,1}) where Dc
  vds = first(visualization_data(trian, ""))
  x = getindex.(get_node_coordinates(vds.grid),1)
  x, zeros(size(x))
end
