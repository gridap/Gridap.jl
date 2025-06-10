
struct VisualizationData
  grid::Grid
  filebase::AbstractString
  celldata
  nodaldata
  function VisualizationData(grid::Grid,filebase::AbstractString;celldata=Dict(),nodaldata=Dict())
    new(grid,filebase,celldata,nodaldata)
  end
end

"""
This function returns an iterable collection (e.g. a Vector) of VisualizationData objects
"""
function visualization_data(args...;kwargs...)
  @abstractmethod
end

struct VisualizationGrid{Dc,Dp} <: Grid{Dc,Dp}
  sub_grid::Grid{Dc,Dp}
  sub_cell_to_cell::AbstractVector{<:Integer}
  cell_to_refpoints::AbstractVector{<:AbstractVector{<:Point}}
end

Geometry.get_reffes(g::VisualizationGrid) = get_reffes(g.sub_grid)
Geometry.get_cell_type(g::VisualizationGrid) = get_cell_type(g.sub_grid)
Geometry.get_node_coordinates(g::VisualizationGrid) = get_node_coordinates(g.sub_grid)
Geometry.get_cell_node_ids(g::VisualizationGrid) = get_cell_node_ids(g.sub_grid)

function visualization_grid(args...;kwargs...)
  @abstractmethod
end

# Triangulations and Grids

function visualization_data(grid::Grid, filebase::AbstractString;celldata=Dict(),nodaldata=Dict())
  (VisualizationData(grid,filebase;celldata=celldata,nodaldata=nodaldata),)
end

function visualization_data(
  trian::Triangulation, filebase::AbstractString;
  order=-1, nsubcells=-1, celldata=Dict(), cellfields=Dict())
  visualization_data(get_grid(trian),trian,filebase;order,nsubcells,celldata,cellfields)
end

function visualization_data(
  ::Grid, trian::Triangulation, filebase::AbstractString;
  order=-1, nsubcells=-1, celldata=Dict(), cellfields=Dict())

  if order == -1 && nsubcells == -1
    # Use the given cells as visualization cells
    f = (reffe) -> UnstructuredGrid(reffe)
  elseif order != -1 && nsubcells == -1
    # Use cells of given order as visualization cells
    f = (reffe) -> UnstructuredGrid(LagrangianRefFE(Float64,get_polytope(reffe),order))
  elseif order == -1 && nsubcells != -1
    # Use linear sub-cells with nsubcells per direction
    f = (reffe) -> UnstructuredGrid(compute_reference_grid(reffe,nsubcells))
  else
    @unreachable "order and nsubcells kw-arguments can not be given at the same time"
  end

  ref_grids = map(f, get_reffes(trian))
  visgrid = visualization_grid(trian,ref_grids)

  cdata = _prepare_cdata(celldata,visgrid.sub_cell_to_cell)
  pdata = _prepare_pdata(trian,cellfields,visgrid.cell_to_refpoints)

  (VisualizationData(visgrid,filebase;celldata=cdata,nodaldata=pdata),)
end

function visualization_grid(trian::Triangulation, ref_grids::AbstractArray{<:UnstructuredGrid})
  cell_to_ctype = collect1d(get_cell_type(trian))
  ctype_to_refpoints = map(get_node_coordinates, ref_grids)
  cell_to_refpoints = CompressedArray(ctype_to_refpoints,cell_to_ctype)
  
  cell_map = get_cell_map(trian)
  cell_to_points = lazy_map(evaluate,cell_map,cell_to_refpoints)
  node_to_coords, cell_to_offset = _prepare_node_to_coords(cell_to_points)

  ctype_to_subcell_to_snodes = map(get_cell_node_ids,ref_grids)
  sub_cell_to_nodes, sub_cell_to_cell = _prepare_sub_cell_to_nodes(
    cell_to_ctype,ctype_to_subcell_to_snodes,cell_to_offset
  )

  ctype_to_reffes = map(get_reffes,ref_grids)
  ctype_to_scell_type = map(get_cell_type,ref_grids)
  sctype_to_reffe, sub_cell_to_sctype = _prepare_sctype_to_reffe(
    cell_to_ctype,ctype_to_reffes,ctype_to_scell_type
  )

  sub_grid = UnstructuredGrid(
    node_to_coords, sub_cell_to_nodes, sctype_to_reffe, sub_cell_to_sctype, NonOriented()
  )

  return VisualizationGrid(sub_grid,sub_cell_to_cell,cell_to_refpoints)
end

# Polytopal Grids

function visualization_data(
  ::Geometry.PolytopalGrid, trian::Triangulation, filebase::AbstractString;
  order=-1, nsubcells=-1, celldata=Dict(), cellfields=Dict())
  @notimplementedif order != -1 "order kw-argument is not implemented for PolytopalGrid"
  @notimplementedif nsubcells != -1 "nsubcells kw-argument is not implemented for PolytopalGrid"

  grid = visualization_grid(get_grid(trian))
  cdata = _prepare_cdata(celldata,Base.OneTo(num_cells(trian)))
  pdata = _prepare_pdata(trian,cellfields,get_cell_coordinates(grid))

  (VisualizationData(grid,filebase;celldata=cdata,nodaldata=pdata),)
end

function visualization_grid(grid::Geometry.PolytopalGrid)
  polytopes = get_polytopes(grid)
  cell_to_coords = map(get_vertex_coordinates,polytopes)
  node_coords = reduce(vcat,cell_to_coords)
  
  ptrs = zeros(Int32,length(polytopes)+1)
  for (i,p) in enumerate(polytopes)
    ptrs[i+1] = num_vertices(p)
  end
  length_to_ptrs!(ptrs)
  data = collect(Int,1:length(node_coords))
  cell_to_nodes = Table(data,ptrs)

  return Geometry.PolytopalGrid(node_coords,cell_to_nodes,polytopes)
end

# Polytope and ReferenceFE

function visualization_data(p::Polytope,filebase)
  map(0:(num_dims(p)-1)) do d
    filebase_d = "$(filebase)_$d"
    grid = Grid(ReferenceFE{d},p)
    VisualizationData(grid,filebase_d)
  end
end

function visualization_data(x::AbstractVector{<:Point}, filebase; kwargs...)
  grid = UnstructuredGrid(x)
  (VisualizationData(grid,filebase; kwargs...),)
end

function visualization_data(reffe::LagrangianRefFE,filebase)
  p = get_polytope(reffe)
  visdata_p = visualization_data(p,filebase)
  node_coords = get_node_coordinates(reffe)
  node_comp_to_dof = get_node_and_comp_to_dof(reffe)
  nodaldata = [
    "dof" => node_comp_to_dof,
    "node" => collect(1:num_nodes(reffe))
  ]
  visdata_n = visualization_data(node_coords, "$(filebase)_nodes"; nodaldata = nodaldata)
  [visdata_p..., visdata_n...]
end

# DiscreteModel

function visualization_data(
  model::DiscreteModel, filebase::AbstractString;labels::FaceLabeling=get_face_labeling(model))
  map(0:num_cell_dims(model)) do d
    grid = Grid(ReferenceFE{d},model)
    cdata = _prepare_cdata(labels,d)
    VisualizationData(grid,"$(filebase)_$d";celldata=cdata)
  end
end

# CellPoint

function visualization_data(x::CellPoint,filename; cellfields=Dict())
  cell_to_points = get_array(x)
  cf(v) = CellField(v,get_triangulation(x),DomainStyle(x))
  nodaldata = Dict((k => evaluate(cf(v),x) for (k,v) in cellfields))
  visualization_data(cell_to_points, filename, nodaldata=nodaldata)
end

function visualization_data(
  cell_to_points::AbstractArray{<:AbstractArray{<:Point}}, filename; 
  celldata=Dict(), nodaldata=Dict()
)
  node_to_point, cell_to_offset = _prepare_node_to_coords(cell_to_points)
  nnodes = length(node_to_point)
  cell_to_nodes = identity_table(Int,Int32,nnodes)
  cell_to_ctype = fill(Int8(1),nnodes)
  ctype_to_reffe = [VERTEX1,]
  node_to_cell = _prepare_node_to_cell(cell_to_offset,nnodes)
  grid = UnstructuredGrid(
    node_to_point,
    cell_to_nodes,
    ctype_to_reffe,
    cell_to_ctype,
    Oriented())
  cdata = _prepare_cdata(celldata,node_to_cell)
  pdata = _prepare_pdata(nodaldata)

  (VisualizationData(grid,filename,celldata=cdata,nodaldata=pdata),)
end

# Private functions

function _prepare_pdata(trian,cellfields,cellpoints)
  x = CellPoint(cellpoints,trian,ReferenceDomain())
  nodaldata = Dict((k => evaluate(CellField(v,trian),x) for (k,v) in cellfields))
  return _prepare_pdata(nodaldata)
end

function _prepare_pdata(nodaldata)
  pdata = Dict()
  for (k,v) in nodaldata
    pdata[k], = _prepare_node_to_coords(v)
  end
  return pdata
end

function _prepare_cdata(celldata,sub_cell_to_cell)
  cdata = Dict()
  for (k,v) in celldata
    cell_to_a = collect(get_array(v))
    T = eltype(cell_to_a)
    sub_cell_to_a = zeros(T,length(sub_cell_to_cell))
    for (sub_cell, cell) in enumerate(sub_cell_to_cell)
      sub_cell_to_a[sub_cell] = cell_to_a[cell]
    end
    cdata[k] = sub_cell_to_a
  end
  @assert ! haskey(cdata,"cell") "cell is a reserved key"
  cdata["cell"] = sub_cell_to_cell
  return cdata
end

function _prepare_cdata(labels::FaceLabeling,d)
  cdata = []
  dface_to_entity = get_face_entity(labels,d)
  for tag in 1:num_tags(labels)
    dface_to_isontag = zeros(Int,num_faces(labels,d))
    for entity in get_tag_entities(labels,tag)
      for (i,e) in enumerate(dface_to_entity)
        if isequal(e,entity)
          dface_to_isontag[i] = tag
        end
      end
    end
    name = get_tag_name(labels,tag)
    push!(cdata, name => dface_to_isontag )
  end
  push!(cdata,"entity" => dface_to_entity)
  return cdata
end

function _prepare_node_to_coords(cell_to_points)
  cell_to_offset = zeros(Int,length(cell_to_points))
  P = eltype(eltype(cell_to_points))
  node_to_coords = P[]

  offset = 0
  cache = array_cache(cell_to_points)
  for cell in eachindex(cell_to_points)
    cell_to_offset[cell] = offset
    points = getindex!(cache,cell_to_points,cell)
    append!(node_to_coords,points)
    offset += length(points)
  end

  return node_to_coords, cell_to_offset
end

function _prepare_sub_cell_to_nodes(
  cell_to_ctype,ctype_to_scell_to_snodes,cell_to_offset)

  n_scells = 0
  for ctype in cell_to_ctype
    n_scells += length(ctype_to_scell_to_snodes[ctype])
  end

  scell = 1
  scell_to_cell = Vector{Int}(undef,n_scells)
  ptrs = Vector{Int32}(undef,n_scells+1)
  for (cell, ctype) in enumerate(cell_to_ctype)
    scell_to_snodes = ctype_to_scell_to_snodes[ctype]
    for snodes in scell_to_snodes
      scell_to_cell[scell] = cell
      ptrs[scell+1] = length(snodes)
      scell += 1
    end
  end
  length_to_ptrs!(ptrs)
  @assert scell == n_scells+1

  data = Vector{Int}(undef,ptrs[end]-1)
  scell = 1
  for (cell, ctype) in enumerate(cell_to_ctype)
    scell_to_snodes = ctype_to_scell_to_snodes[ctype]
    offset = cell_to_offset[cell]
    for snodes in scell_to_snodes
      for snode in snodes
        data[ptrs[scell]] = snode + offset
        ptrs[scell] += 1
      end
      scell += 1
    end
  end
  rewind_ptrs!(ptrs)
  @assert scell == n_scells+1

  scell_to_nodes = Table(data,ptrs)
  return scell_to_nodes, scell_to_cell
end

function _prepare_sctype_to_reffe(cell_to_ctype,ctype_to_r_to_reffe,ctype_to_scell_to_r)

  i_to_reffe = vcat(ctype_to_r_to_reffe...)
  u_to_reffe, i_to_u = _find_unique_with_indices(i_to_reffe)

  n_ctypes = length(ctype_to_r_to_reffe)
  ctype_to_nreffes = map(length,ctype_to_r_to_reffe)

  n_scells = 0 
  for ctype in cell_to_ctype
    n_scells += length(ctype_to_scell_to_r[ctype])
  end
  
  i = 1
  ctype_to_r_to_u = Vector{Vector{Int}}(undef,n_ctypes)
  for ctype in eachindex(ctype_to_r_to_reffe)
    nreffes = ctype_to_nreffes[ctype]
    r_to_u = zeros(Int,nreffes)
    for r in 1:nreffes
      r_to_u[r] = i_to_u[i]
      i += 1
    end
    ctype_to_r_to_u[ctype] = r_to_u
  end

  scell = 1
  scell_to_u = Vector{Int8}(undef,n_scells)
  for ctype in cell_to_ctype
    scell_to_r = ctype_to_scell_to_r[ctype]
    r_to_u = ctype_to_r_to_u[ctype]
    for r in scell_to_r
      scell_to_u[scell] = r_to_u[r]
      scell += 1
    end
  end

  return u_to_reffe, scell_to_u
end

function _prepare_node_to_cell(cell_to_offset,nnodes)
  node_to_cell = zeros(Int,nnodes)
  ncells = length(cell_to_offset)
  for cell in 1:ncells
    nini = cell_to_offset[cell]+1
    if cell < ncells
      nend = cell_to_offset[cell+1]
    else
      nend = nnodes
    end
    for node in nini:nend
      node_to_cell[node] = cell
    end
  end
  node_to_cell
end
