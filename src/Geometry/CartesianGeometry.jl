module CartesianGeometry

# Dependencies of this module

using StaticArrays: SVector, MVector, @SVector
using Gridap
using Gridap.Helpers

# Functionality provided

export CartesianGrid
export CartesianDiscreteModel
import Base: size, getindex, IndexStyle
import Gridap: points, cells, celltypes, cellorders
import Gridap: Grid, FaceLabels
import Gridap: UnstructuredGrid
import Gridap: FlexibleUnstructuredGrid
import Gridap: GridGraph
import Gridap: FullGridGraph

struct CartesianGrid{D} <: Grid{D,D}
  dim_to_limits::NTuple{D,NTuple{2,Float64}}
  dim_to_ncells::NTuple{D,Int}
  extrusion::NTuple{D,Int}
  order:: Int
  map
end

const identity = (x -> x)

function CartesianGrid(;partition::NTuple{D,Int},domain=nothing,
                        order::Int=1,map=identity) where D
  _cartesiangrid(partition,domain,order,map)
end

function points(self::CartesianGrid)
  dim_to_npoint = tuple([ i+1 for i in self.dim_to_ncells ]...)
  CartesianGridPoints(self.map,self.dim_to_limits,dim_to_npoint)
end

cells(self::CartesianGrid) = CartesianGridCells(self.dim_to_ncells)

celltypes(self::CartesianGrid) = ConstantCellValue(self.extrusion,prod(self.dim_to_ncells))

cellorders(self::CartesianGrid) = ConstantCellValue(self.order,prod(self.dim_to_ncells))

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
struct CartesianDiscreteModel{D} <: DiscreteModel{D}
  cgrid::CartesianGrid{D}
  ugrid::UnstructuredGrid
  gridgraph::FullGridGraph
end

function CartesianDiscreteModel(; args...)
  cgrid = CartesianGrid(; args...)
  ugrid = UnstructuredGrid(cgrid)
  gridgraph = FullGridGraph(ugrid)
  CartesianDiscreteModel(cgrid, ugrid, gridgraph)
end

Grid(model::CartesianDiscreteModel{D},::Val{D}) where D = model.cgrid

function Grid(model::CartesianDiscreteModel{D},::Val{Z}) where {D,Z}

  ugrid = model.ugrid
  gridgraph = model.gridgraph
  face_to_vertices = connections(gridgraph,Z,0)
  nfaces = length(face_to_vertices)
  fcode = tuple([HEX_AXIS for i in 1:Z]...)

  order = 1
  @notimplementedif order != (cellorders(model.cgrid)).value

  _points = points(ugrid)
  _cells_data = face_to_vertices.data
  _cells_ptrs = face_to_vertices.ptrs
  _ctypes = ConstantCellValue(fcode,nfaces)
  _corders = ConstantCellValue(order,nfaces)

  UnstructuredGrid(
    _points, _cells_data, _cells_ptrs, _ctypes, _corders)
end

GridGraph(model::CartesianDiscreteModel) = model.gridgraph

FullGridGraph(model::CartesianDiscreteModel) = model.gridgraph

#@fverdugo precompute this result
function FaceLabels(model::CartesianDiscreteModel{D}) where D
  dim_to_face_to_geolabel = Vector{Vector{Int}}(undef,D+1)
  dim_to_offset = _generate_dim_to_offset(D)
  interior_id = dim_to_offset[end]+1
  boundary_id = -1
  for d in 0:(D-1)
    face_to_cells = connections(model.gridgraph,d,D)
    cell_to_faces = connections(model.gridgraph,D,d)
    offset = dim_to_offset[d+1]
    face_to_geolabel = _generate_pre_geolabel(
      face_to_cells,
      cell_to_faces,
      offset,
      interior_id,boundary_id,d,D)
    dim_to_face_to_geolabel[d+1] = face_to_geolabel
  end
  for d in 0:(D-2)
    for j in (d+1):(D-1)
      dface_to_jfaces = connections(model.gridgraph,d,j)
      dface_to_geolabel = dim_to_face_to_geolabel[d+1]
      jface_to_geolabel = dim_to_face_to_geolabel[j+1]
      _fix_dface_geolabels!(
        dface_to_geolabel,
        jface_to_geolabel,
        dface_to_jfaces,
        interior_id,boundary_id)
    end
  end
  _ncells = ncells(model.cgrid)
  dim_to_face_to_geolabel[end] = ConstantCellValue(interior_id,_ncells)
  phys_labels = [ [i] for i in 1:interior_id ]
  push!(phys_labels,[i for i in 1:(interior_id-1)])
  tag_to_name = ["physical_tag_$i" for i in 1:interior_id-1]
  push!(tag_to_name,"interior")
  push!(tag_to_name,"boundary")
  FaceLabels(dim_to_face_to_geolabel, phys_labels, tag_to_name)
end

# Helpers

function _generate_dim_to_offset(D)
  code = tuple(fill(HEX_AXIS,D)...)
  polytope = Polytope(code)
  dim_to_offset = zeros(Int,D+1)
  for d in 1:D
    dim_to_offset[d+1] = dim_to_offset[d] + num_nfaces(polytope,d-1)
  end
  dim_to_offset
end

function _generate_pre_geolabel(
  face_to_cells,
  cell_to_faces,
  offset,
  interior_id,
  boundary_id,d,D)

  nfaces = length(face_to_cells)
  face_to_geolabel = fill(interior_id,nfaces)

  max_ncells_around = 2^(D-d)

  _generate_pre_geolabel_kernel!(
    face_to_geolabel,
    face_to_cells.data,
    face_to_cells.ptrs,
    cell_to_faces.data,
    cell_to_faces.ptrs,
    offset,boundary_id,max_ncells_around)

  face_to_geolabel
end

function _generate_pre_geolabel_kernel!(
  face_to_geolabel,
  face_to_cells_data,
  face_to_cells_ptrs,
  cell_to_faces_data,
  cell_to_faces_ptrs,
  offset,boundary_id,max_ncells_around)

  nfaces = length(face_to_geolabel)
  for face in 1:nfaces
    a = face_to_cells_ptrs[face]-1
    ncells_around = face_to_cells_ptrs[face+1] - (a+1)
    if ncells_around == 1
      icell_around = 1
      cell = face_to_cells_data[a+icell_around]
      b = cell_to_faces_ptrs[cell]-1
      nlfaces = cell_to_faces_ptrs[cell+1] - (b+1)
      for lface in 1:nlfaces
        face2 = cell_to_faces_data[b+lface]
        if face == face2
          face_to_geolabel[face] = lface + offset
          break
        end
      end
    elseif ncells_around != max_ncells_around
      face_to_geolabel[face] = boundary_id
    end
  end

end

function _fix_dface_geolabels!(
  dface_to_geolabel,
  jface_to_geolabel,
  dface_to_jfaces,
  interior_id,boundary_id)

  for (dface, jfaces) in enumerate(dface_to_jfaces)
    if dface_to_geolabel[dface] != boundary_id
      continue
    end
    for jface in jfaces
      geolabel = jface_to_geolabel[jface]
      if geolabel != interior_id && geolabel != boundary_id
        dface_to_geolabel[dface] = geolabel
        break
      end
    end
  end

end

function _cartesiangrid(partition::NTuple{D,Int},domain,order,map) where D
  if domain === nothing
    _domain = [ i*j for i in ones(D) for j in 0:1 ]
    # _domain = [ i*(-1)^j for i in ones(D) for j in 1:2 ]
  else
    _domain = domain
  end
  dim_to_limits = tuple([(_domain[2*i-1],_domain[2*i]) for i in 1:D ]...)
  extrusion = tuple(fill(HEX_AXIS,D)...)
  dim_to_ncells = partition
  @notimplementedif order != 1
  CartesianGrid{D}(dim_to_limits,dim_to_ncells,extrusion,order,map)
end

struct CartesianGridPoints{D} <: IndexCellValue{Point{D,Float64},D}
  map
  dim_to_limits::NTuple{D,NTuple{2,Float64}}
  dim_to_npoint::NTuple{D,Int}
end

size(self::CartesianGridPoints) = self.dim_to_npoint

IndexStyle(::Type{CartesianGridPoints{D}} where D) = IndexCartesian()

function getindex(self::CartesianGridPoints{D}, I::Vararg{Int, D}) where D
  p = zero(MVector{D,Float64})
  @inbounds for d in 1:D
    xa = self.dim_to_limits[d][1]
    xb = self.dim_to_limits[d][2]
    p[d] =  xa + (I[d]-1)*(xb-xa)/(self.dim_to_npoint[d]-1)
  end
  # map(Point(p)
  # self.map != identity ? self.map(Point(p)) : Point(p)
  self.map(Point(p))
end

struct CartesianGridCells{D,L} <: IndexCellArray{Int,1,SVector{L,Int},D}
  dim_to_ncell::SVector{D,Int}
end

function CartesianGridCells(dim_to_ncell::NTuple{D,Int}) where D
  CartesianGridCells{D,2^D}(dim_to_ncell)
end

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
  ps = Array{Point{D,Float64},1}(undef,(length(points(grid)),))
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

end # module
