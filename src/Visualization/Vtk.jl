
#"""
#"""
#function writevtk(
#  trian::Grid, filebase; celldata=Dict(), nodaldata=Dict())
#  write_vtk_file(trian,filebase,celldata=celldata,nodaldata=nodaldata)
#end

"""
    writevtk(reffe::LagrangianRefFE,filebase)
"""
function writevtk(reffe::LagrangianRefFE,filebase)

  p = get_polytope(reffe)
  writevtk(p,filebase)

  node_coords = get_node_coordinates(reffe)
  node_comp_to_dof = get_node_and_comp_to_dof(reffe)
  nodaldata = [
    "dof" => node_comp_to_dof,
    "node" => collect(1:num_nodes(reffe))]
  writevtk(node_coords, "$(filebase)_nodes"; nodaldata = nodaldata)

end

"""
    writevtk(x::AbstractVector{<:Point}, filebase; kwargs...)
"""
function writevtk(x::AbstractVector{<:Point}, filebase; kwargs...)
  vtkfile = createvtk(x,filebase; kwargs...)
  outfiles = vtk_save(vtkfile)
end

"""
    createvtk(x::AbstractVector{<:Point}, filebase; kwargs...)
"""
function createvtk(x::AbstractVector{<:Point}, filebase; kwargs...)
  grid = UnstructuredGrid(x)
  create_vtk_file(grid,filebase; kwargs...)
end

"""
    writevtk(p::Polytope,filebase)
"""
function writevtk(p::Polytope,filebase)
  for d in 0:(num_dims(p)-1)
    grid = Grid(ReferenceFE{d},p)
    write_vtk_file(grid,"$(filebase)_$d")
  end
end

"""
    writevtk(model::DiscreteModel,filebase)
    writevtk(model::DiscreteModel, labels::FaceLabeling, filebase)
"""
function writevtk(model::DiscreteModel,filebase)
  labels = get_face_labeling(model)
  writevtk(model,labels,filebase)
end

function writevtk(model::DiscreteModel, labels::FaceLabeling, filebase)
  for d in 0:num_cell_dims(model)
    grid = Grid(ReferenceFE{d},model)
    cdat = _prepare_cdata_model(labels,d)
    write_vtk_file(grid,"$(filebase)_$d";celldata=cdat)
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

"""

    write_vtk_file(
      trian::Grid,
      filebase;
      celldata=Dict(),
      nodaldata=Dict())

Low level entry point to vtk. Other vtk-related routines in Gridap eventually call this one.

"""
function write_vtk_file(
  trian::Grid, filebase; celldata=Dict(), nodaldata=Dict())
  vtkfile = create_vtk_file(trian, filebase, celldata=celldata, nodaldata=nodaldata)
  outfiles = vtk_save(vtkfile)
end

"""

    create_vtk_file(
      trian::Grid,
      filebase;
      celldata=Dict(),
      nodaldata=Dict())

Low level entry point to vtk. Other vtk-related routines in Gridap eventually call this one.
This function only creates the vtkFile, without writing to disk.

"""
function create_vtk_file(
  trian::Grid, filebase; celldata=Dict(), nodaldata=Dict())

  points = _vtkpoints(trian)
  cells = _vtkcells(trian)
  vtkfile = vtk_grid(filebase, points, cells, compress=false)

  for (k,v) in celldata
    vtk_cell_data(vtkfile, _prepare_data(v), k)
  end
  for (k,v) in nodaldata
    vtk_point_data(vtkfile, _prepare_data(v), k)
  end

  return vtkfile
end

function _vtkpoints(trian)
  D = num_point_dims(trian)
  x = get_node_coordinates(trian)
  xflat = collect(x)
  reshape(reinterpret(Float64,xflat),(D,length(x)))
end

function _vtkcells(trian)

  type_to_reffe = get_reffes(trian)
  cell_to_type = get_cell_type(trian)
  type_to_vtkid = map(get_vtkid, type_to_reffe)
  type_to_vtknodes = map(get_vtknodes, type_to_reffe)
  cell_to_nodes = get_cell_nodes(trian)
  cache = array_cache(cell_to_nodes)

  _generate_vtk_cells(
    cache,
    cell_to_nodes,
    cell_to_type,
    type_to_vtkid,
    type_to_vtknodes)
end

function _generate_vtk_cells(
  cache,
  cell_to_nodes,
  cell_to_type,
  type_to_vtkid,
  type_to_vtknodes)

  V = eltype(cell_to_nodes)
  meshcells = MeshCell{WriteVTK.VTKCellTypes.VTKCellType,V}[]

  d = _vtkcelltypedict()

  cells = 1:length(cell_to_type)
  for cell in cells

    t = cell_to_type[cell]
    vtkid = type_to_vtkid[t]
    vtknodes = type_to_vtknodes[t]

    nodes = getindex!(cache,cell_to_nodes,cell)
    meshcell = MeshCell(d[vtkid], nodes[vtknodes])

    push!(meshcells,meshcell)

  end

  meshcells

end

_prepare_data(v) = v

function _prepare_data(v::AbstractArray{<:MultiValue})
  a = collect(v)
  reinterpret(a)
end

function _prepare_data(v::AbstractArray{<:VectorValue{2,T}}) where T
  a = collect(v)
  b = reshape(reinterpret(T,a),(2,length(a)))
  z = zeros((1,size(b,2)))
  vcat(b,z)
end

"""
"""
function get_vtkid(reffe::LagrangianRefFE)
  basis = get_prebasis(reffe)
  p = get_polytope(reffe)
  get_vtkid(p,basis)
end

"""
"""
function get_vtknodes(reffe::LagrangianRefFE)
  basis = get_prebasis(reffe)
  p = get_polytope(reffe)
  get_vtknodes(p,basis)
end

function get_vtkid(p::ExtrusionPolytope, basis::MonomialBasis)
  exponents = get_exponents(basis)
  vtkid, _ = _vtkinfo_extrusion_polytope(p,exponents)
  vtkid
end

function get_vtknodes(p::ExtrusionPolytope, basis::MonomialBasis)
  exponents = get_exponents(basis)
  _, vtknodes = _vtkinfo_extrusion_polytope(p,exponents)
  vtknodes
end

function get_vtkid(p::SerendipityPolytope,basis::MonomialBasis)
  get_vtkid(p.hex,basis)
end

function get_vtknodes(p::SerendipityPolytope,basis::MonomialBasis)
  get_vtknodes(p.hex,basis)
end

function _vtkinfo_extrusion_polytope(p,exponents)

  # Taken from the vtk specification
  # https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

  n_nodes = length(exponents)

  if p == VERTEX
    if n_nodes == 1
      vtkid = 1
      vtknodes = [1,]
    else
      @notimplemented
    end

  elseif p == SEGMENT
    if n_nodes == 2
      vtkid = 3
      vtknodes = [1,2]
    elseif n_nodes == 3
      vtkid = 21
      vtknodes = [1,2,3]
    else
      @notimplemented
    end

  elseif p == TRI
    if n_nodes == 3
      vtkid = 5
      vtknodes = [1,2,3]
    elseif n_nodes == 6
      vtkid = 22
      vtknodes = [1,2,3,4,6,5]
    elseif n_nodes == 10
      vtkid = 69
      vtknodes = [1,2,3,4,5,8,9,7,6,10]
    elseif n_nodes == 15
      vtkid = 69
      vtknodes = [1,2,3,4,5,6,10,11,12,9,8,7,13,14,15]
    elseif n_nodes == 21
      vtkid = 69
      vtknodes = [1,2,3,4,5,6,7,12,13,14,15,11,10,9,8,16,18,21,17,20,19]
    else
      @notimplemented
    end

  elseif p == QUAD
    if n_nodes == 4
      vtkid = 9
      vtknodes = [1,2,4,3]
    elseif n_nodes == 8
      vtkid = 23
      vtknodes = [1,2,4,3,5,8,6,7]
    elseif n_nodes == 9
      vtkid = 28
      vtknodes = [1,2,4,3,5,8,6,7,9]
    elseif n_nodes == 16
      vtkid = 70
      vtknodes = [1,2,4,3,5,6,11,12,7,8,9,10,13,14,15,16]
    elseif n_nodes == 25
      vtkid = 70
      vtknodes = [1,2,4,3,5,6,7,14,15,16,8,9,10,11,12,13,17,18,19,20,21,22,23,24,25]
    else
      @notimplemented
    end

  elseif p == TET
    if n_nodes == 4
      vtkid = 10
      vtknodes = [1,2,3,4]
    elseif n_nodes == 10
      vtkid = 24
                 #0,1,2,3,4,5,6,7,8,9
      vtknodes = [1,2,3,4,5,7,6,8,9,10]
    elseif n_nodes == 20
      vtkid = 71
      vtknodes = [1,2,3,4,5,6,9,10,7,8,11,12,13,14,15,16,18,20,19,17]
    else
      @notimplemented
    end

  elseif p == HEX
    if n_nodes == 8
      vtkid = 12
      vtknodes = [1,2,4,3,5,6,8,7]
    elseif n_nodes == 27
      vtkid = 72
      vtknodes = [1,2,4,3,5,6,8,7,9,14,10,13,11,16,12,15,17,18,19,20,25,26,23,24,21,22,27]
    elseif n_nodes == 64
      vtkid = 72
      vtknodes = [1,2,4,3,5,6,8,7,9,10,19,20,11,12,17,18,13,14,23,24,15,16,21,22,25,26,27,
                  28,29,30,31,32,49,50,51,52,53,54,55,56,41,42,43,44,45,46,47,48,33,34,35,
                  36,37,38,39,40,57,58,59,60,61,62,63,64]
    else
      @notimplemented
    end

  elseif p == WEDGE
    if n_nodes == 6
      vtkid = 13
      vtknodes = [1,3,2,4,6,5]
    else
      @notimplemented
    end

  elseif p == PYRAMID
    if n_nodes == 5
      vtkid = 14
      vtknodes = [1,2,4,3,5]
    else
      @notimplemented
    end

  else
    @notimplemented "vtkid not implemented for given ExtrusionPolytope"
  end

  (vtkid, vtknodes)
end

function _vtkcelltypedict()
  d = Dict{Int,WriteVTK.VTKCellTypes.VTKCellType}()
  d[VTK_VERTEX.vtk_id] = VTK_VERTEX
  d[VTK_LINE.vtk_id] = VTK_LINE
  d[VTK_LINE.vtk_id] = VTK_LINE
  d[VTK_TRIANGLE.vtk_id] = VTK_TRIANGLE
  d[VTK_QUAD.vtk_id] = VTK_QUAD
  d[VTK_TETRA.vtk_id] = VTK_TETRA
  d[VTK_HEXAHEDRON.vtk_id] = VTK_HEXAHEDRON
  d[VTK_WEDGE.vtk_id] = VTK_WEDGE
  d[VTK_QUADRATIC_QUAD.vtk_id] = VTK_QUADRATIC_QUAD
  d[VTK_BIQUADRATIC_QUAD.vtk_id] = VTK_BIQUADRATIC_QUAD
  d[VTK_QUADRATIC_TRIANGLE.vtk_id] = VTK_QUADRATIC_TRIANGLE
  d[VTK_QUADRATIC_TETRA.vtk_id] = VTK_QUADRATIC_TETRA
  d[VTK_QUADRATIC_EDGE.vtk_id] = VTK_QUADRATIC_EDGE
  d[VTK_QUADRATIC_HEXAHEDRON.vtk_id] = VTK_QUADRATIC_HEXAHEDRON
  d[VTK_PYRAMID.vtk_id] = VTK_PYRAMID
  d[VTK_LAGRANGE_TRIANGLE.vtk_id] = VTK_LAGRANGE_TRIANGLE
  d[VTK_LAGRANGE_QUADRILATERAL.vtk_id] = VTK_LAGRANGE_QUADRILATERAL
  d[VTK_LAGRANGE_TETRAHEDRON.vtk_id] = VTK_LAGRANGE_TETRAHEDRON
  d[VTK_LAGRANGE_HEXAHEDRON.vtk_id] = VTK_LAGRANGE_HEXAHEDRON
  #d[VTK_BIQUADRATIC_HEXAHEDRON.vtk_id] = VTK_BIQUADRATIC_HEXAHEDRON
  d
end

"""
"""
function writevtk(trian::Triangulation, filebase; order=-1, nsubcells=-1, celldata=Dict(), cellfields=Dict())
  vtkfile = createvtk(trian,filebase, order=order, nsubcells=nsubcells, celldata=celldata, cellfields=cellfields)
  outfiles = vtk_save(vtkfile)
end

"""
"""
function createvtk(trian::Triangulation, filebase; order=-1, nsubcells=-1, celldata=Dict(), cellfields=Dict())
  visdata = visualization_data(trian; order=order,
    nsubcells=nsubcells, celldata=celldata, cellfields=cellfields)
  create_vtk_file(visdata.grid,filebase,celldata=visdata.celldata,nodaldata=visdata.nodaldata)
end

function _prepare_pdata(trian,cellfields,samplingpoints)
  pdata = Dict()
  for (k,v) in cellfields
    _v = CellField(v,trian)
    pdata[k], = _prepare_node_to_coords(evaluate(_v,samplingpoints))
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

"""
    writevtk(
      cell_to_points::AbstractArray{<:AbstractArray{<:Point}},
      filename;
      celldata=Dict(),
      nodaldata=Dict())
"""
function writevtk(
  cell_to_points::AbstractArray{<:AbstractArray{<:Point}},
  filename; celldata=Dict(), nodaldata=Dict())
  vtkfile = createvtk(cell_to_points, filename, celldata=celldata, nodaldata=nodaldata)
  outfiles = vtk_save(vtkfile)
end


"""
    createvtk(
      cell_to_points::AbstractArray{<:AbstractArray{<:Point}},
      filename;
      celldata=Dict(),
      nodaldata=Dict())
"""
function createvtk(
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
    Val{true}())

  cdata = _prepare_cdata(celldata,node_to_cell)
  pdata = _prepare_pdata_for_cell_points(nodaldata)

  create_vtk_file(grid,filename,celldata=cdata,nodaldata=pdata)

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
