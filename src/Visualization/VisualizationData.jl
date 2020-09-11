struct VisualizationGrid{Dc,Dp,T} <: Grid{Dc,Dp}
  sub_grid::UnstructuredGrid{Dc,Dp,T,false}
  sub_cell_to_cell::Vector{Int}
  cell_to_refpoints::CompressedArray{Vector{Point{Dc,T}},1,Vector{Vector{Point{Dc,T}}},Vector{Int8}}
end

get_reffes(g::VisualizationGrid) = get_reffes(g.sub_grid)

get_cell_type(g::VisualizationGrid) = get_cell_type(g.sub_grid)

get_node_coordinates(g::VisualizationGrid) = get_node_coordinates(g.sub_grid)

get_cell_nodes(g::VisualizationGrid) = get_cell_nodes(g.sub_grid)

function VisualizationGrid(trian::Triangulation, ref_grids::Vector{<:UnstructuredGrid})

  cell_to_ctype = collect1d(get_cell_type(trian))
  ctype_to_refpoints = map(get_node_coordinates, ref_grids)
  cell_to_refpoints = CompressedArray(ctype_to_refpoints,cell_to_ctype)
  cell_map = get_cell_map(trian)
  cell_to_points = evaluate(cell_map, cell_to_refpoints)

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
    Val{false}())

  VisualizationGrid(sub_grid,sub_cell_to_cell,cell_to_refpoints)

end

struct VisualizationData
  grid::VisualizationGrid
  celldata
  nodaldata
end

"""
"""
function visualization_data(trian::Triangulation; order=-1, nsubcells=-1, celldata=Dict(), cellfields=Dict())

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

  return VisualizationData(visgrid, cdata, pdata)
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
  ϕ = get_cell_map(trian)
  q = GenericCellPoint(samplingpoints)
  x = ϕ(q)
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
    cell_to_a = collect(v)
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
