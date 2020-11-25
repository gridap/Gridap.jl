
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

# Visualizing Triangulation, cellfields, etc.

function visualization_data(
  trian::Triangulation, filebase::AbstractString;
  order=-1, nsubcells=-1, celldata=Dict(), cellfields=Dict())

  if order == -1 && nsubcells == -1
    # Use the given cells as visualization cells
    f = (reffe) -> UnstructuredGrid(reffe)
  elseif order != -1 && nsubcells == -1
    # Use cells of given order as visualization cells
    f = (reffe) -> UnstructuredGrid(LagrangianRefFE(Float64,get_polytope(reffe),order))
  elseif order == -1 && nsubcells != -1
    # Use use linear sub-cells with nsubcells per direction
    f = (reffe) -> UnstructuredGrid(compute_reference_grid(reffe,nsubcells))
  else
    @unreachable "order and nsubcells kw-arguments can not be given at the same time"
  end

  ref_grids = map(f, get_reffes(trian))
  visgrid = VisualizationGrid(trian,ref_grids)

  cdata = _prepare_cdata(celldata,visgrid.sub_cell_to_cell)
  pdata = _prepare_pdata(trian,cellfields,visgrid.cell_to_refpoints)

  (VisualizationData(visgrid,filebase;celldata=cdata,nodaldata=pdata),)
end

function _prepare_node_to_coords(cell_to_points)
  cell_to_offset = zeros(Int,length(cell_to_points))
  P = eltype(eltype(cell_to_points))
  node_to_coords = P[]
  cache = array_cache(cell_to_points)
  _prepare_node_to_coords_fill!(node_to_coords,cell_to_offset,cache,cell_to_points)
  (node_to_coords, cell_to_offset)
end

function _prepare_node_to_coords_fill!(node_to_coords,cell_to_offset,cache,cell_to_points)
  offset = 0
  for cell in 1:length(cell_to_points)
    cell_to_offset[cell] = offset
    points = getindex!(cache,cell_to_points,cell)
    for point in points
      push!(node_to_coords, point)
      offset += 1
    end
  end
end

function _prepare_sub_cell_to_nodes(
  cell_to_ctype,ctype_to_scell_to_snodes,cell_to_offset)

  sub_cell_to_nodes = Table(Int[],Int32[1])
  sub_cell_to_cell = Int[]

  for (cell, ctype) in enumerate(cell_to_ctype)
    scell_to_snodes = ctype_to_scell_to_snodes[ctype]
    offset = cell_to_offset[cell]
    for (scell, snodes) in enumerate(scell_to_snodes)
      push!(sub_cell_to_cell,cell)
      push!(sub_cell_to_nodes.ptrs,length(snodes))
      for snode in snodes
        node = snode+offset
        push!(sub_cell_to_nodes.data,node)
      end
    end
  end

  length_to_ptrs!(sub_cell_to_nodes.ptrs)

  (sub_cell_to_nodes, sub_cell_to_cell)
end

function _prepare_sctype_to_reffe(ctype_to_r_to_reffe,ctype_to_scell_to_r,cell_to_ctype)

  i_to_reffe = vcat(ctype_to_r_to_reffe...)
  u_to_reffe, i_to_u = _find_unique_with_indices(i_to_reffe)
  ctype_to_nreffes = map(length,ctype_to_r_to_reffe)
  sub_cell_to_u = _prepare_sub_cell_to_u(
    ctype_to_r_to_reffe,ctype_to_nreffes,ctype_to_scell_to_r,cell_to_ctype,i_to_u)

  (u_to_reffe, sub_cell_to_u)
end

function _prepare_sub_cell_to_u(
  ctype_to_r_to_reffe,ctype_to_nreffes,ctype_to_scell_to_r,cell_to_ctype,i_to_u)

  i = 1
  ctype_to_r_to_u = Vector{Int}[]
  for (ctype, r_to_reffes) in enumerate(ctype_to_r_to_reffe)
    nreffes = ctype_to_nreffes[ctype]
    r_to_u = zeros(Int,nreffes)
    for r in 1:nreffes
      r_to_u[r] = i_to_u[i]
      i +=1
    end
    push!(ctype_to_r_to_u,r_to_u)
  end

  sub_cell_to_u = Int8[]
  for ctype in cell_to_ctype
    scell_to_r = ctype_to_scell_to_r[ctype]
    r_to_u = ctype_to_r_to_u[ctype]
    for r in scell_to_r
      push!(sub_cell_to_u,r_to_u[r])
    end
  end

  sub_cell_to_u

end

function _prepare_pdata(trian,cellfields,samplingpoints)
  x = CellPoint(samplingpoints,trian,ReferenceDomain())
  pdata = Dict()
  for (k,v) in cellfields
    _v = CellField(v,trian)
    pdata[k], = _prepare_node_to_coords(evaluate(_v,x))
  end
  pdata
end

function _prepare_cdata(celldata,sub_cell_to_cell)
  cdata = Dict()
  for (k,v) in celldata
    cell_to_a = collect(get_array(v))
    T = eltype(cell_to_a)
    sub_cell_to_a = zeros(T,length(sub_cell_to_cell))
    _fill_sub_cell_to_a!(sub_cell_to_a,cell_to_a,sub_cell_to_cell)
    cdata[k] = sub_cell_to_a
  end
  k2 = "cell"
  @assert ! haskey(cdata,k2) "cell is a reserved key"
  cdata[k2] = sub_cell_to_cell
  cdata
end

function _fill_sub_cell_to_a!(sub_cell_to_a,cell_to_a,sub_cell_to_cell)
  for (sub_cell, cell) in enumerate(sub_cell_to_cell)
    sub_cell_to_a[sub_cell] = cell_to_a[cell]
  end
end

struct VisualizationGrid{Dc,Dp} <: Grid{Dc,Dp}
  sub_grid::UnstructuredGrid{Dc,Dp}
  sub_cell_to_cell::AbstractVector{<:Integer}
  cell_to_refpoints::AbstractVector{<:AbstractVector{<:Point}}
end

get_reffes(g::VisualizationGrid) = get_reffes(g.sub_grid)

get_cell_type(g::VisualizationGrid) = get_cell_type(g.sub_grid)

get_node_coordinates(g::VisualizationGrid) = get_node_coordinates(g.sub_grid)

get_cell_nodes(g::VisualizationGrid) = get_cell_nodes(g.sub_grid)

function VisualizationGrid(trian::Triangulation, ref_grids::AbstractArray{<:UnstructuredGrid})

  cell_to_ctype = collect1d(get_cell_type(trian))
  ctype_to_refpoints = map(get_node_coordinates, ref_grids)
  cell_to_refpoints = CompressedArray(ctype_to_refpoints,cell_to_ctype)
  cell_map = get_cell_map(trian)
  cell_to_points = lazy_map(evaluate,cell_map, cell_to_refpoints)

  node_to_coords, cell_to_offset = _prepare_node_to_coords(cell_to_points)

  ctype_to_scell_to_snodes = map(get_cell_nodes,ref_grids)

  sub_cell_to_nodes, sub_cell_to_cell = _prepare_sub_cell_to_nodes(
    cell_to_ctype,ctype_to_scell_to_snodes,cell_to_offset)

  ctype_to_reffes = map(get_reffes,ref_grids)
  ctype_to_scell_type = map(get_cell_type,ref_grids)
  sctype_to_reffe, sub_cell_to_sctype = _prepare_sctype_to_reffe(
    ctype_to_reffes,ctype_to_scell_type,cell_to_ctype)

  sub_grid = UnstructuredGrid(
    node_to_coords,
    sub_cell_to_nodes,
    sctype_to_reffe,
    sub_cell_to_sctype,
    NonOriented())

  VisualizationGrid(sub_grid,sub_cell_to_cell,cell_to_refpoints)

end

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
    "node" => collect(1:num_nodes(reffe))]
  visdata_n = visualization_data(node_coords, "$(filebase)_nodes"; nodaldata = nodaldata)
  [visdata_p..., visdata_n...]
end

function visualization_data(
  model::DiscreteModel, filebase::AbstractString;labels::FaceLabeling=get_face_labeling(model))
  map(0:num_cell_dims(model)) do d
    grid = Grid(ReferenceFE{d},model)
    cdat = _prepare_cdata_model(labels,d)
    VisualizationData(grid,"$(filebase)_$d";celldata=cdat)
  end
end

function _prepare_cdata_model(labels,d)
  dface_to_entity = get_face_entity(labels,d)
  cdat = []
  for tag in 1:num_tags(labels)
    dface_to_isontag = zeros(Int,num_faces(labels,d))
    for entity in get_tag_entities(labels,tag)
      _set_entity!(dface_to_isontag,dface_to_entity,entity,tag)
    end
    name = get_tag_name(labels,tag)
    push!(cdat, name => dface_to_isontag )
  end
  push!(cdat,"entity" => dface_to_entity)
  cdat
end

function _set_entity!(dface_to_isontag,dface_to_entity,entity,tag)
  for i in 1:length(dface_to_entity)
    if dface_to_entity[i] == entity
      dface_to_isontag[i] = tag
    end
  end
end

function visualization_data(
  cell_to_points::AbstractArray{<:AbstractArray{<:Point}},
  filename; celldata=Dict(), nodaldata=Dict())

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
  pdata = _prepare_pdata_for_cell_points(nodaldata)

  (VisualizationData(grid,filename,celldata=cdata,nodaldata=pdata),)
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

function _prepare_pdata_for_cell_points(nodaldata)
  pdata = Dict()
  for (k,v) in nodaldata
    pdata[k], = _prepare_node_to_coords(v)
  end
  pdata
end

function visualization_data(x::CellPoint,filename; cellfields=Dict())
  nodaldata = Dict((
    k=>evaluate(CellField(v,get_triangulation(x),DomainStyle(x)),x) for (k,v) in cellfields ))
  cell_to_points = get_array(x)
  visualization_data(cell_to_points, filename, nodaldata=nodaldata)
end

