
"""
"""
function writevtk(
  trian::Grid, filebase; celldata=Dict(), nodaldata=Dict())
  write_vtk_file(trian,filebase,celldata=celldata,nodaldata=nodaldata)
end

"""
    writevtk(reffe::NodalReferenceFE,filebase)
"""
function writevtk(reffe::NodalReferenceFE,filebase)

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
  grid = UnstructuredGrid(x)
  write_vtk_file(grid,filebase; kwargs...)
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

  points = _vtkpoints(trian)
  cells = _vtkcells(trian)
  vtkfile = vtk_grid(filebase, points, cells, compress=false)

  for (k,v) in celldata
    vtk_cell_data(vtkfile, _prepare_data(v), k)
  end
  for (k,v) in nodaldata
    vtk_point_data(vtkfile, _prepare_data(v), k)
  end

  outfiles = vtk_save(vtkfile)
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
  meshcells = MeshCell{V}[]

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

function _prepare_data(v::AbstractArray{<:VectorValue{D,T}}) where {D,T}
  a = collect(v)
  reshape(reinterpret(T,a),(D,length(a)))
end

function _prepare_data(v::AbstractArray{<:VectorValue{2,T}}) where T
  a = collect(v)
  b = reshape(reinterpret(T,a),(2,length(a)))
  z = zeros((1,size(b,2)))
  vcat(b,z)
end

function _prepare_data(v::AbstractArray{<:TensorValue{D,T}}) where {D,T}
  a = collect(v)
  reshape(reinterpret(T,a),(D*D,length(a)))
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
    else
      @notimplemented
    end

  elseif p == HEX
    if n_nodes == 8
      vtkid = 12
      vtknodes = [1,2,4,3,5,6,8,7]
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
  #d[VTK_BIQUADRATIC_HEXAHEDRON.vtk_id] = VTK_BIQUADRATIC_HEXAHEDRON
  d
end

