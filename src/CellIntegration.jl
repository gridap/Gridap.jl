module CellIntegration

export IntegrationMesh
export geomap, cellcoordinates, cellbasis, ncells
export integrate

using Numa.Helpers
using Numa.FieldValues
using Numa.Polynomials
using Numa.CellValues
using Numa.CellFunctions
using Numa.CellQuadratures
using Numa.Quadratures

import Numa: evaluate, gradient
import Numa: cellfield
import Numa.Geometry: celltypes

"""
Minimal interface for a mesh used for numerical integration
"""
abstract type IntegrationMesh{Z,D} end

function cellcoordinates(::IntegrationMesh{Z,D})::CellPoints{D} where {Z,D}
 @abstractmethod
end

function cellbasis(::IntegrationMesh{Z,D})::CellBasis{Z,Float64} where {Z,D}
  @abstractmethod
end

"""
Returns the tuple uniquely identifying the Polytope of each cell
"""
function celltypes(::IntegrationMesh{Z,D})::CellValue{NTuple{Z}} where {Z,D}
  @abstractmethod
end

function geomap(self::IntegrationMesh)
  coords = cellcoordinates(self)
  basis = cellbasis(self)
  expand(basis,coords)
end

function ncells(self::IntegrationMesh)
  coords = cellcoordinates(self)
  length(coords)
end

function integrate(cellfun::CellFunction{Point{D},1,T,N},phi::CellGeomap{D,Z},quad::CellQuadrature{D}) where {D,Z,T,N}
  z = coordinates(quad)
  w = weights(quad)
  f = evaluate(cellfun,z)
  j = evaluate(gradient(phi),z)
  cellsum( f*(meas(j)*w), dim=N )
end

function integrate(cellfun::CellFunction{Point{D},1},mesh::IntegrationMesh{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = geomap(mesh)
  integrate(cellfun,phi,quad)
end

function integrate(fun::Function,mesh::IntegrationMesh{D,Z},quad::CellQuadrature{D}) where {D,Z}
  phi = geomap(mesh)
  cellfun = compose(fun,phi)
  integrate(cellfun,phi,quad)
end

function cellfield(mesh::IntegrationMesh,fun::Function)
  phi = geomap(mesh)
  compose(fun,phi)
end

function cellfield(mesh::IntegrationMesh{D,Z},fun::Function,u::CellField{Z}) where {D,Z}
  phi = geomap(mesh)
  compose(fun,phi,u)
end

# @fverdugo Avoid code repetition with a generated function

using Numa: flatten
using Numa.Polytopes
using Numa.Geometry

function refgrid end

function refgrid(::Val{(HEX_AXIS,)},nref::Int)
  n = 2^nref
  CartesianGrid(domain=(-1.0,1.0),partition=(n,))
end

function refgrid(::Val{(HEX_AXIS,HEX_AXIS)},nref::Int)
  n = 2^nref
  CartesianGrid(domain=(-1.0,1.0,-1.0,1.0),partition=(n,n))
end

function refgrid(::Val{(HEX_AXIS,HEX_AXIS,HEX_AXIS)},nref::Int)
  n = 2^nref
  CartesianGrid(domain=(-1.0,1.0,-1.0,1.0,-1.0,1.0),partition=(n,n,n))
end

# @fverdugo Move this to another place if we move the definition
# of IntegrationMesh
function visgrid(self::IntegrationMesh;nref=0)
  grid, coarsecells, samplingpoints = _prepare_grid(celltypes(self),geomap(self),nref)
  VisualizationGrid(grid,coarsecells,samplingpoints)
end

function _prepare_grid(celltypes::CellValue{NTuple{Z,Int}},phi::CellGeomap{Z,D},nref) where {Z,D}
  @notimplemented
end

function _prepare_grid(ctypes::ConstantCellValue{NTuple{Z,Int}},phi::CellGeomap{Z,D},nref) where {Z,D}
  refgrid = _prepare_refgrid(ctypes,nref)
  samplingpoints = _prepare_samplingpoints(refgrid,ctypes)
  ps, offsets = _prepare_points(samplingpoints,points(refgrid),phi)
  data, ptrs, coarsecells = _prepare_cells(refgrid,offsets)
  ts = _prepare_celltypes(length(ctypes),celltypes(refgrid))
  grid = UnstructuredGrid(ps,data,ptrs,ts)
  (grid, coarsecells, samplingpoints)
end

function _prepare_refgrid(ctypes,nref)
  extrusion = celldata(ctypes)
  refgrid(Val(extrusion),nref)
end

function _prepare_samplingpoints(refgrid,ctypes)
  refpoints = flatten(collect(points(refgrid)))
  ConstantCellArray(refpoints,length(ctypes))
end

function _prepare_points(samplingpoints,refpoints,phi::CellGeomap{Z,D}) where {Z,D}
  xe = evaluate(phi,samplingpoints)
  offsets = Array{Int,1}(undef,(length(xe),))
  ps = Array{Point{D},1}(undef,(length(xe)*length(refpoints)))
  _fill_ps_and_offsets!(ps,offsets,xe)
  (ps, offsets)
end

function _fill_ps_and_offsets!(ps,offsets,xe)
  k = 1
  for (n,x) in enumerate(xe)
    offsets[n] = k - 1
    for xi in x
      @inbounds ps[k] = xi
      k += 1
    end
  end
end

function _prepare_cells(refgrid,offsets)
  refcells = cells(refgrid)
  ptrs = Array{Int,1}(undef,(1+length(offsets)*length(refcells),))
  coarsecells = Array{Int,1}(undef,(length(offsets)*length(refcells),))
  _fill_ptrs_and_coarsecells!(ptrs,coarsecells,refcells,length(offsets))
  length_to_ptrs!(ptrs)
  data = Array{Int,1}(undef,(ptrs[end]-1,))
  _fill_data!(data,offsets,refcells)
  (data, ptrs, CellValueFromArray(coarsecells))
end

function _fill_ptrs_and_coarsecells!(ptrs,coarsecells,refcells,ncells)
  k = 1
  for cell in 1:ncells
    for refnodes in refcells
      @inbounds coarsecells[k] = cell
      k += 1
      @inbounds ptrs[k] = length(refnodes)
    end
  end
end

function _fill_data!(data,offsets,refcells)
  k = 1
  for offset in offsets
    for refnodes in refcells
      for node in refnodes
        @inbounds data[k] = node + offset
        k += 1
      end
    end
  end
end

function _prepare_celltypes(ncells,refcelltypes::ConstantCellValue)
  refextrusion = celldata(refcelltypes)
  fill(refextrusion,ncells*length(refcelltypes) )
end

end # module CellIntegration
