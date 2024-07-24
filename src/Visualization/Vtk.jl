"""
"""
function writevtk(args...;compress=false,append=true,ascii=false,vtkversion=:default,kwargs...)
  map(visualization_data(args...;kwargs...)) do visdata
    write_vtk_file(
      visdata.grid,visdata.filebase,celldata=visdata.celldata,nodaldata=visdata.nodaldata,
      compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
    )
  end
end

"""
"""
function createvtk(args...;compress=false,append=true,ascii=false,vtkversion=:default,kwargs...)
  v = visualization_data(args...;kwargs...)
  @notimplementedif length(v) != 1
  visdata = first(v)
  create_vtk_file(
    visdata.grid,visdata.filebase,celldata=visdata.celldata,nodaldata=visdata.nodaldata,
    compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
  )
end

"""
"""
function createpvd(args...;kwargs...)
  paraview_collection(args...;kwargs...)
end

function createpvd(parts::Nothing,args...;kwargs...)
  createpvd(args...;kwargs...)
end

function createpvd(f,parts::Nothing,args...;kwargs...)
  createpvd(f,args...;kwargs...)
end

"""
"""
function savepvd(pvd::WriteVTK.CollectionFile)
  vtk_save(pvd)
end

"""

    write_vtk_file(
      trian::Grid,
      filebase;
      celldata=Dict(),
      nodaldata=Dict(),
      vtk_kwargs...
      )

Low level entry point to vtk. Other vtk-related routines in Gridap eventually call this one.
The optional WriteVTK kwargs `vtk_kwargs` are passed to the `vtk_grid` constructor.
"""
function write_vtk_file(
  trian::Grid, filebase; celldata=Dict(), nodaldata=Dict(),
  compress=false, append=true, ascii=false, vtkversion=:default
)
  vtkfile = create_vtk_file(
    trian, filebase, celldata=celldata, nodaldata=nodaldata,
    compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
  )
  outfiles = vtk_save(vtkfile)
end

"""

    create_vtk_file(
      trian::Grid,
      filebase;
      celldata=Dict(),
      nodaldata=Dict(),
      vtk_kwargs...
    )

Low level entry point to vtk. Other vtk-related routines in Gridap eventually call this one.
This function only creates the vtkFile, without writing to disk.

The optional WriteVTK kwargs `vtk_kwargs` are passed to the `vtk_grid` constructor.
"""
function create_vtk_file(
  trian::Grid, filebase; celldata=Dict(), nodaldata=Dict(), 
  compress=false, append=true, ascii=false, vtkversion=:default
)

  points = _vtkpoints(trian)
  cells = _vtkcells(trian)
  vtkfile = vtk_grid(
    filebase, points, cells,
    compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
  )

  if num_cells(trian)>0
    for (k,v) in celldata
      vtk_cell_data(vtkfile, _prepare_data(v), k)
    end
    for (k,v) in nodaldata
      vtk_point_data(vtkfile, _prepare_data(v), k)
    end
  end

  return vtkfile
end

function create_pvtk_file(
  trian::Grid, filebase; part, nparts, ismain=(part==1), celldata=Dict(), nodaldata=Dict(),
  compress=false, append=true, ascii=false, vtkversion=:default
)

  points = _vtkpoints(trian)
  cells = _vtkcells(trian)
  vtkfile = pvtk_grid(
    filebase, points, cells;part=part, nparts=nparts, ismain=ismain,
    compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
  )

  if num_cells(trian) > 0
    for (k, v) in celldata
      vtkfile[k, VTKCellData()] = _prepare_data(v)
    end
    for (k, v) in nodaldata
      vtkfile[k, VTKPointData()] = _prepare_data(v)
    end
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
  cell_to_nodes = get_cell_node_ids(trian)
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
