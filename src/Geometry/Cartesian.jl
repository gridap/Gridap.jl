module Cartesian

# Dependencies of this module

using StaticArrays: SVector, MVector, @SVector
using UnstructuredGrids: UGrid
using UnstructuredGrids: generate_dual_connections
using UnstructuredGrids: generate_cell_to_faces
using UnstructuredGrids: connections
using Numa.Helpers
using Numa.FieldValues
using Numa.Polytopes
using Numa.Meshes
using Numa.Geometry
using Numa.Geometry.Unstructured
using Numa.CellValues

# Functionality provided 

export CartesianGrid
import Base: size, getindex, IndexStyle
import Numa.CellValues: cellsize
import Numa.Geometry: points, cells, celltypes, cellorders, gridgraph
import Numa.Geometry: Grid, GridGraph, NFaceLabels, boundarylabels
import Numa.Geometry.Unstructured: UnstructuredGrid
import Numa.Geometry.Unstructured: FlexibleUnstructuredGrid
export CartesianDiscreteModel

struct CartesianGrid{D} <: Grid{D,D}
  dim_to_limits::NTuple{D,NTuple{2,Float64}}
  dim_to_ncells::NTuple{D,Int}
  extrusion::NTuple{D,Int}
  order:: Int
end

function CartesianGrid(;partition::NTuple{D,Int},domain=nothing,order::Int=1) where D
  _cartesiangrid(partition,domain,order)
end

function points(self::CartesianGrid)
  dim_to_npoint = tuple([ i+1 for i in self.dim_to_ncells ]...)
  CartesianGridPoints(self.dim_to_limits,dim_to_npoint)
end

cells(self::CartesianGrid) = CartesianGridCells(self.dim_to_ncells)

celltypes(self::CartesianGrid) = ConstantCellValue(self.extrusion,prod(self.dim_to_ncells))

cellorders(self::CartesianGrid) = ConstantCellValue(self.order,prod(self.dim_to_ncells))

function gridgraph(self::CartesianGrid)
  #fverdugo this is a temporary implementation
  nparts = [i for i in self.dim_to_ncells]
  mesh = StructHexMesh(nparts)
  GridGraphFromData(mesh.cellvefs,mesh.vefcells)
end

"""
Construct an `UnstructuredGrid` from a `CartesianGrid`
"""
function UnstructuredGrid(grid::CartesianGrid{D}) where D
  ps = _compute_points(grid)
  ts = celltypes(grid)
  os = cellorders(grid)
  data, ptrs = _compute_cells(grid,UnstructuredGrid{D,D})
  UnstructuredGrid(ps,data,ptrs,ts,os)
end

"""
Construct an `FlexibleUnstructuredGrid` from a `CartesianGrid`
"""
function FlexibleUnstructuredGrid(grid::CartesianGrid{D}) where D
  ps = _compute_points(grid)
  ts = celltypes(grid)
  os = cellorders(grid)
  cs = _compute_cells(grid,FlexibleUnstructuredGrid{D,D})
  FlexibleUnstructuredGrid(ps,cs,ts,os)
end

"""
DiscreteModel associated with a CartesianGrid
"""
struct CartesianDiscreteModel{
  D,U<:UGrid,G<:NewGridGraph{D}} <: DiscreteModel{D}

  cgrid::CartesianGrid{D}
  ugrid::U
  gridgraph::G
end

function CartesianDiscreteModel(; args...)
  cgrid = CartesianGrid(; args...)
  grid = UnstructuredGrid(cgrid)
  ugrid = UGrid(grid)
  gridgraph = _setup_grid_graph(ugrid,celldim(cgrid))
  CartesianDiscreteModel(cgrid, ugrid, gridgraph)
end

Grid(model::CartesianDiscreteModel{D},::Val{D}) where D = model.cgrid

function Grid(model::CartesianDiscreteModel{D},::Val{Z}) where {D,Z}
  @notimplemented
end

GridGraph(model::CartesianDiscreteModel{D},::Val{D}) where D = model.gridgraph

function GridGraph(::CartesianDiscreteModel{D},::Val{Z}) where {D,Z}
  @notimplemented
end

function NFaceLabels(::CartesianDiscreteModel{D}) where D
  @notimplemented
end

function boundarylabels(::CartesianDiscreteModel)
  @notimplemented
end

# Helpers

function _setup_grid_graph(ugrid,D)

  cell_to_vertices = connections(ugrid)
  vertex_to_cells = generate_dual_connections(cell_to_vertices)
  dim_to_cell_to_vefs = Vector{IndexCellArray{Int,1}}(undef,D)
  dim_to_vef_to_cells = Vector{IndexCellArray{Int,1}}(undef,D)

  d=0
  dim_to_cell_to_vefs[d+1] = CellVectorFromDataAndPtrs(
    cell_to_vertices.list,cell_to_vertices.ptrs)

  dim_to_vef_to_cells[d+1] = CellVectorFromDataAndPtrs(
    vertex_to_cells.list,vertex_to_cells.ptrs)

  for d=1:(D-1)

    cell_to_dfaces = generate_cell_to_faces(d,ugrid,vertex_to_cells)
    dface_to_cells = generate_dual_connections(cell_to_dfaces)

    dim_to_cell_to_vefs[d+1] = CellVectorFromDataAndPtrs(
      cell_to_dfaces.list,cell_to_dfaces.ptrs)

    dim_to_vef_to_cells[d+1] = CellVectorFromDataAndPtrs(
    dface_to_cells.list,dface_to_cells.ptrs)

  end

  NewGridGraph(dim_to_cell_to_vefs, dim_to_vef_to_cells)

end

function _cartesiangrid(partition::NTuple{D,Int},domain,order) where D
  if domain === nothing
    _domain = [ i*(-1)^j for i in ones(D) for j in 1:2 ]
  else
    _domain = domain
  end
  dim_to_limits = tuple([(_domain[2*i-1],_domain[2*i]) for i in 1:D ]...)
  extrusion = tuple(fill(HEX_AXIS,D)...)
  dim_to_ncells = partition
  @notimplementedif order != 1
  CartesianGrid{D}(dim_to_limits,dim_to_ncells,extrusion,order)
end

struct CartesianGridPoints{D} <: IndexCellValue{Point{D},D}
  dim_to_limits::NTuple{D,NTuple{2,Float64}}
  dim_to_npoint::NTuple{D,Int}
end

size(self::CartesianGridPoints) = self.dim_to_npoint

IndexStyle(::Type{CartesianGridPoints{D}} where D) = IndexCartesian()

function getindex(self::CartesianGridPoints{D}, I::Vararg{Int, D}) where D
  p = zero(MPoint{D})
  @inbounds for d in 1:D
    xa = self.dim_to_limits[d][1]
    xb = self.dim_to_limits[d][2]
    p[d] =  xa + (I[d]-1)*(xb-xa)/(self.dim_to_npoint[d]-1)
  end
  Point{D}(p)
end

struct CartesianGridCells{D,L} <: IndexCellArray{Int,1,SVector{L,Int},D}
  dim_to_ncell::SVector{D,Int}
end

function CartesianGridCells(dim_to_ncell::NTuple{D,Int}) where D
  CartesianGridCells{D,2^D}(dim_to_ncell)
end

cellsize(self::CartesianGridCells{D,L}) where {D,L} = (L,)

size(self::CartesianGridCells) = self.dim_to_ncell.data

IndexStyle(::Type{CartesianGridCells{D,L}} where {D,L}) = IndexCartesian()

function getindex(self::CartesianGridCells{D,L}, I::Vararg{Int, D}) where {D,L}
  dim_to_ngpoint = 1 .+ self.dim_to_ncell
  dim_to_nlpoint = @SVector fill(2,D)
  offset = @SVector fill(1,D)
  pointgids = LinearIndices(dim_to_ngpoint.data)
  cellpointlids = CartesianIndices(dim_to_nlpoint.data)
  cellgid = CartesianIndex(I...) - CartesianIndex(offset.data)
  cellpointgids = cellpointlids .+ cellgid
  ids = zero(MVector{L,Int})
  @inbounds for (l,pgid) in enumerate(cellpointgids)
    ids[l] = pointgids[pgid]
  end
  SVector{L,Int}(ids)
end

function _compute_points(grid::CartesianGrid{D}) where D
  ps = Array{Point{D},1}(undef,(length(points(grid)),))
  for (i,xi) in enumerate(points(grid))
    ps[i] = xi
  end
  ps
end

function _compute_cells(grid::CartesianGrid{D},::Type{FlexibleUnstructuredGrid{D,D}}) where D
  cs = [ Array{Int,1}(undef,(2^D,)) for i in 1:length(cells(grid)) ]
  for (i,ci) in enumerate(cells(grid))
    cs[i] .= ci
  end
  cs
end

function _compute_cells(grid::CartesianGrid{D},::Type{UnstructuredGrid{D,D}}) where D
  ptrs = fill(2^D,(length(cells(grid))+1,))
  length_to_ptrs!(ptrs)
  data = zeros(Int,ptrs[end]-1)
  _fill_cell_data!(data,cells(grid))
  (data, ptrs)
end

function _fill_cell_data!(data,cells)
  k = 1
  for v in cells
    for vi in v
      @inbounds data[k] = vi
      k +=1
    end
  end
end

end # module Cartesian
