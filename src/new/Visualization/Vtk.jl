

"""

Low level entry point to vtk. Other, vtk-related routines in Gridap eventually call this one.

The reference FEs in the triangulation need to implement the following queries

- [`get_vtkid(reffe::NodalReferenceFE)`](@ref)
- [`get_vtknodes(reffe::NodalReferenceFE)`](@ref)

A default implementation is available for `LagrangianRefFE`. It is based on the following
queries on the underlying polytope

- [`get_vtkid(p::Polytope,basis::MonomialBasis)`](@ref)
- [`get_vtknodes(p::Polytope,basis::MonomialBasis)`](@ref)

They are implemented for `ExtrusionPolytope`.

"""
function writevtk_file(
  trian::ConformingTriangulation, filebase; celldata=Dict(), nodaldata=Dict())

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

function _prepare_data(v::AbstractArray{<:VectorValue{D}}) where D
  a = collect(v)
  reshape(reinterpret(Float64,a),(D,length(a)))
end

function _prepare_data(v::AbstractArray{<:VectorValue{2}})
  a = collect(v)
  b = reshape(reinterpret(Float64,a),(2,length(a)))
  z = zeros((1,size(b,2)))
  vcat(b,z)
end

function _prepare_data(v::AbstractArray{<:TensorValue{D}}) where D
  a = collect(v)
  reshape(reinterpret(Float64,a),(D*D,length(a)))
end

"""
"""
function get_vtkid(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
"""
function get_vtknodes(reffe::NodalReferenceFE)
  @abstractmethod
end

function get_vtkid(reffe::LagrangianRefFE)
  basis = get_prebasis(reffe)
  p = get_polytope(reffe)
  get_vtkid(p,basis)
end

function get_vtknodes(reffe::LagrangianRefFE)
  basis = get_prebasis(reffe)
  p = get_polytope(reffe)
  get_vtknodes(p,basis)
end

"""
    get_vtkid(p::Polytope,basis::MonomialBasis) -> Int

Given a polytope `p` and a monomial basis, returns an integer with its vtk identifier.
Overloading of this function is needed only in order to visualize the underlying polytope
with Paraview.
"""
function get_vtkid(p::Polytope,basis::MonomialBasis)
  @abstractmethod
end

"""
    get_vtknodes(p::Polytope,basis::MonomialBasis) -> Vector{Int}

Given a polytope `p` and monomial basis, returns a vector of integers representing a permutation of the
polytope vertices required to relabel the vertices according the criterion adopted in
Paraview.
Overloading of this function is needed only in order to visualize the underlying polytope
with Paraview.
"""
function get_vtknodes(p::Polytope,basis::MonomialBasis)
  @abstractmethod
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
    else
      @notimplemented
    end

  elseif p == TRI
    if n_nodes == 3
      vtkid = 5
      vtknodes = [1,2,3]
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
  d[VTK_QUADRATIC_QUAD.vtk_id] = VTK_QUADRATIC_QUAD
  d[VTK_BIQUADRATIC_QUAD.vtk_id] = VTK_BIQUADRATIC_QUAD
  d
end

# Visualization of a LagrangianRefFE



