module Draft

using Gridap.Helpers
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Arrays
using Gridap.Visualization
using FillArrays

export CellBoundary

struct CellBoundaryGrid{Dc,Dp,P,A,B,C,D} <: Grid{Dc,Dp}
  parent::P
  node_coord::A
  cell_lface_nodes::B
  ftype_freffe::C
  cell_lface_ftype::D
  function CellBoundaryGrid(grid::Grid)
    D = num_cell_dims(grid)
    Dp = num_point_dims(grid)
    node_coord = get_node_coordinates(grid)
    ctype_reffe = get_reffes(grid)
    @notimplementedif length(ctype_reffe) != 1
    reffe = first(ctype_reffe)
    freffes = get_reffaces(ReferenceFE{D-1},reffe)
    @notimplementedif length(freffes) != 1
    ftype_freffe = [first(freffes),]
    ctype_lface_ftype = map(reffe->get_face_type(reffe,D-1),ctype_reffe)
    ctype_lface_to_lnodes = map(reffe->get_face_nodes(reffe,D-1),ctype_reffe)
    cell_ctype = get_cell_type(grid)
    cell_nodes = get_cell_node_ids(grid)
    cell_lface_ftype = expand_cell_data(ctype_lface_ftype,cell_ctype)
    cell_lface_nodes = lazy_map(cell_ctype,cell_nodes) do ctype,lnode_node
      # This can be heavily optimized
      lface_to_lnodes = ctype_lface_to_lnodes[ctype]
      lface_nodes = [ lnode_node[lnodes] for lnodes in  lface_to_lnodes]
      lface_nodes
    end
    A = typeof(node_coord)
    B = typeof(cell_lface_nodes)
    C = typeof(ftype_freffe)
    E = typeof(ctype_lface_ftype)
    P = typeof(grid)
    new{D-1,Dp,P,A,B,C,E}(
      grid,node_coord,cell_lface_nodes,ftype_freffe,cell_lface_ftype)
  end
end

Geometry.get_node_coordinates(a::CellBoundaryGrid) = a.node_coord
Geometry.get_cell_node_ids(a::CellBoundaryGrid) = a.cell_lface_nodes
Geometry.get_reffes(a::CellBoundaryGrid) = a.ftype_freffe
Geometry.get_cell_type(a::CellBoundaryGrid) = a.cell_lface_ftype

# The following ones are a bit "hacky"

function Geometry.get_cell_coordinates(trian::CellBoundaryGrid)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_node_ids(trian)
  lazy_map(Broadcasting(Broadcasting(Reindex(node_to_coords))),cell_to_nodes)
end

function Geometry.get_cell_map(trian::CellBoundaryGrid)
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_shapefuns = get_cell_shapefuns(trian)
  lazy_map(Broadcasting(linear_combination),cell_to_coords,cell_to_shapefuns)
end

struct CellBoundaryTriangulation{Dc,Dp,A,B} <: Triangulation{Dc,Dp}
  model::A
  grid::B
  function CellBoundaryTriangulation(model::DiscreteModel)
    A = typeof(model)
    D = num_cell_dims(model)
    grid = CellBoundaryGrid(get_grid(model))
    B = typeof(grid)
    new{D-1,D,A,B}(model,grid)
  end
end

CellBoundary(args...) = CellBoundaryTriangulation(args...)

Geometry.get_background_model(a::CellBoundaryTriangulation) = a.model
Geometry.get_grid(a::CellBoundaryTriangulation) = a.grid

struct CellBoundaryGlue{A,B}
  tcell_lface_mface::A
  tcell_lface_mface_map::B
end

function Geometry.is_change_possible(sglue::FaceToFaceGlue,tglue::CellBoundaryGlue)
  true
end

function Geometry.get_glue(trian::CellBoundaryTriangulation{D},::Val{D}) where D
  model = get_background_model(trian)
  topo = get_grid_topology(model)
  cell_lface_face = get_faces(topo,D+1,D)
  pgrid = trian.grid.parent
  ctype_creffe = get_reffes(pgrid)
  ctype_lface_map = map(ctype_creffe) do reffe
    poly = get_polytope(reffe)
    fill(GenericField(identity),num_faces(poly,D))
  end
  cell_ctype = get_cell_type(pgrid)
  cell_lface_map = expand_cell_data(ctype_lface_map,cell_ctype)
  CellBoundaryGlue(cell_lface_face,cell_lface_map)
end

function Geometry.get_glue(trian::CellBoundaryTriangulation{d},::Val{D}) where {d,D}
  if d+1 != D
    return nothing
  end
  pgrid = trian.grid.parent
  ctype_reffe = get_reffes(pgrid)
  cell_ctype = get_cell_type(pgrid)
  ncells = length(cell_ctype)
  cell_cell = IdentityVector(ncells)
  ctype_nlfaces = map(ctype_reffe) do reffe
    poly = get_polytope(reffe)
    num_faces(poly,d)
  end
  # Avoid allocations here
  tcell_lface_mface = lazy_map(cell_ctype,cell_cell) do ctype, cell
    nlfaces = ctype_nlfaces[ctype]
    fill(cell,nlfaces)
  end
  tcell_lface_mface_map = nothing # TODO: This is the most difficult part
  CellBoundaryGlue(tcell_lface_mface,tcell_lface_mface_map)
end

function CellData.change_domain_ref_ref(
  a::CellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::CellBoundaryGlue)
  @notimplemented # TODO (very important)
end

function CellData.change_domain_phys_phys(
  a::CellField,ttrian::Triangulation,sglue::FaceToFaceGlue,tglue::CellBoundaryGlue)
  sface_to_field = get_data(a)
  mface_to_sface = sglue.mface_to_tface
  tcell_lface_mface = tglue.tcell_lface_mface
  mface_to_field = extend(sface_to_field,mface_to_sface)
  # TODO this can be optimized
  tface_to_field = lazy_map(tcell_lface_mface) do lface_mface
    mface_to_field[lface_mface]
  end
  CellData.similar_cell_field(a,tface_to_field,ttrian,PhysicalDomain())
end

# TODO this needs to be optimized (very important)
function Arrays.evaluate!(
  cache,
  f::AbstractVector,
  x::AbstractVector{<:AbstractVector{<:Point}})
  @check length(f) == length(x)
  evaluate.(f,x)
end

function Visualization.visualization_data(
  a::CellBoundaryTriangulation,
  filebase::AbstractString;
  offset=0,
  cellfields=Dict())

  model = get_background_model(a)
  grid = get_grid(a)
  D = num_cell_dims(model)
  node_coord = get_node_coordinates(grid)
  cell_lface_lfnode_node = get_cell_node_ids(grid)
  ncells = length(cell_lface_lfnode_node)
  parent = grid.parent
  cell_nodes = get_cell_node_ids(parent)

  P = eltype(node_coord)
  fnode_coord = P[]
  face_fnodes = Vector{Int32}[]
  fnode = Int32(0)
  for cell in 1:ncells
    lnode_node = cell_nodes[cell] # Allocation here
    lnode_coord = node_coord[lnode_node]
    Xm = sum(lnode_coord) / length(lnode_coord)
    lface_lfnode_node = cell_lface_lfnode_node[cell] # Allocation here
    for lfnode_node in lface_lfnode_node
      fnodes = Int32[] # Allocation here
      for node in lfnode_node
         Xf = node_coord[node]
         coord = Xf + offset*(Xm-Xf)
         push!(fnode_coord,coord) # Allocation here
         fnode += Int32(1)
         push!(fnodes,fnode) # Allocation here
      end
      push!(face_fnodes,fnodes) # Allocation here
    end
  end

  function compute_pdata(f)
    x = get_cell_points(a)
    cell_lface_node_val = f(x)
    T = eltype(eltype(eltype(cell_lface_node_val)))
    vals = zeros(T,length(fnode_coord))
    i = 0
    # To be optimized
    for lface_node_val in cell_lface_node_val
      for node_val in lface_node_val
        for val in node_val
          i += 1
          vals[i] = val
        end
      end
    end
    vals
  end

  pdata = Dict()
  for (k,v) in cellfields
    pdata[k] = compute_pdata(v)
  end

  ftype_freffe = get_reffes(grid)
  @notimplementedif length(ftype_freffe) != 1
  freffe = first(ftype_freffe)
  freffes = [freffe,]
  nfaces = length(face_fnodes)
  face_ftype = ones(Int8,nfaces)

  fgrid = UnstructuredGrid(
    fnode_coord,Table(face_fnodes),freffes,face_ftype)
  (VisualizationData(fgrid,filebase;nodaldata=pdata),)
end

end # module

