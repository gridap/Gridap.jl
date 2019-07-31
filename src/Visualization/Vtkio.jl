module Vtkio

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery

using WriteVTK
using WriteVTK.VTKCellTypes: VTK_VERTEX
using WriteVTK.VTKCellTypes: VTK_LINE
using WriteVTK.VTKCellTypes: VTK_TRIANGLE
using WriteVTK.VTKCellTypes: VTK_QUAD
using WriteVTK.VTKCellTypes: VTK_TETRA
using WriteVTK.VTKCellTypes: VTK_HEXAHEDRON

# Functionality given by this module

export writevtk
import Gridap: cells, points, celltypes
import UnstructuredGrids: writevtk

"""
Write a Grid object into a vtk file
"""
function writevtk(grid::Grid,filebase;celldata=Dict(),pointdata=Dict())
  _writevtk(grid,filebase,celldata,pointdata)
end

"""
Write a CellPoints object into vtk file
"""
function writevtk(points::CellPoints,filebase;celldata=Dict(),pointdata=Dict())
  _writevtk(points,filebase,celldata,pointdata)
end

"""
Write a CellValue{Point{D}} object into vtk file
"""
function writevtk(points::CellValue{Point{D,Float64}} where D,filebase;celldata=Dict(),pointdata=Dict())
  _writevtk(points,filebase,celldata,pointdata)
end

"""
Write an Triangulation object into vtk
"""
function writevtk(trian::Triangulation,filebase;nref=0,celldata=Dict(),cellfields=Dict())
  _writevtk(trian,filebase,nref,celldata,cellfields)
end

"""
Write a Polytope object into vtk
"""
function writevtk(polytope::Polytope,filebase)
  refcell = RefCell(polytope)
  writevtk(refcell,filebase)
end

"""
Write a DiscreteModel into vtk
"""
function writevtk(model::DiscreteModel{D},filebase) where D
  labels = FaceLabels(model)
  for d in 0:D
    f = "$(filebase)_$(d)"
    grid = Grid(model,d)
    cdat = _prepare_cdata_model(labels,d)
    writevtk(grid,f,celldata=cdat)
  end
end

# Helpers

function _prepare_cdata_model(labels,d)
  dface_to_label = labels_on_dim(labels,d)
  ndfaces = length(dface_to_label)
  cdat = []
  for tag in 1:ntags(labels)
    dface_to_isontag = zeros(Int,ndfaces)
    for label in labels_on_tag(labels,tag)
      _set_label!(dface_to_isontag,dface_to_label,label,tag)
    end
    name = name_from_tag(labels,tag)
    push!(cdat, name => dface_to_isontag )
  end
  push!(cdat,"geolabel" => dface_to_label)
  cdat
end

function _set_label!(dface_to_isontag,dface_to_label,label,tag)
  for i in 1:length(dface_to_label)
    if dface_to_label[i] == label
      dface_to_isontag[i] = tag
    end
  end
end

function _writevtk(grid::Grid,filebase,celldata,pointdata)
  points = _vtkpoints(grid)
  cells = _vtkcells(grid)
  vtkfile = vtk_grid(filebase, points, cells, compress=false)
  for (k,v) in celldata
    vtk_cell_data(vtkfile, _prepare_data(v), k)
  end
  for (k,v) in pointdata
    vtk_point_data(vtkfile, _prepare_data(v), k)
  end
  outfiles = vtk_save(vtkfile)
end

function _vtkpoints(grid::Grid{D}) where D
  x = points(grid)
  xflat = collect(x)
  reshape(reinterpret(Float64,xflat),(D,length(x)))
end

# @fverdugo this allocates a lot of small objects
# Not very crucial since it is for visualization
# but it would be nice to have a better way
function _vtkcells(grid::Grid)
  types = _vtkcelltypedict()
  nodes = _vtkcellnodesdict()
  c = celltypes(grid)
  n = cells(grid)
  _check_order(cellorders(grid))
  [ MeshCell(types[encode_extrusion(ci)], ni[nodes[encode_extrusion(ci)]])
     for (ci,ni) in zip(c,n) ]
end

_check_order(co) = @notimplemented

function _check_order(co::ConstantCellValue{Int})
  o = co.value
  @notimplementedif o != 1
end

"""
Encodes the tuple defining a Polytope into an integer
"""
function encode_extrusion(extrusion::NTuple{Z,Int}) where Z
  k = 0
  for (i,v) in enumerate(extrusion)
    k += v*3^i
  end
  k
end

"""
Decodes an integer into a tuple defining a Polytope
"""
function decode_extrusion(i::Int,::Val{Z}) where Z
  @notimplemented
end


# @fverdugo This can also be done by dispatching on value
"""
Generates the lookup table (as a Dict) in order to convert between
Gridap Polytope identifiers into VTK cell type identifiers
"""
function _vtkcelltypedict()
  d = Dict{Int,WriteVTK.VTKCellTypes.VTKCellType}()
  h = HEX_AXIS
  t = TET_AXIS
  d[encode_extrusion(())] = VTK_VERTEX
  d[encode_extrusion((t,))] = VTK_LINE
  d[encode_extrusion((h,))] = VTK_LINE
  d[encode_extrusion((t,t))] = VTK_TRIANGLE
  d[encode_extrusion((h,h))] = VTK_QUAD
  d[encode_extrusion((t,t,t))] = VTK_TETRA
  d[encode_extrusion((h,h,h))] = VTK_HEXAHEDRON
  d
end

# @fverdugo This can also be done by dispatching on value
"""
Generates the lookup table (as a Dict) in order to convert between
Gridap Polytope corner numbering into VTK corner numbering
"""
function _vtkcellnodesdict()
  d = Dict{Int,Vector{Int}}()
  h = HEX_AXIS
  t = TET_AXIS
  d[encode_extrusion(())] = [1,]
  d[encode_extrusion((t,))] = [1,2]
  d[encode_extrusion((h,))] = [1,2]
  d[encode_extrusion((t,t))] = [1,2,3]
  d[encode_extrusion((h,h))] = [1,2,4,3]
  d[encode_extrusion((t,t,t))] = [1,2,3,4]
  d[encode_extrusion((h,h,h))] = [1,2,4,3,5,6,8,7]
  d
end

function _writevtk(points::CellPoints,filebase,celldata,pointdata)
  grid, p_to_cell = _cellpoints_to_grid(points)
  pdat = _prepare_pointdata(pointdata)
  k = "cellid"
  @assert ! haskey(pdat,k)
  pdat[k] = p_to_cell
  writevtk(grid,filebase,pointdata=pdat)
end

function _writevtk(points::CellValue{Point{D,Float64}} where D,filebase,celldata,pointdata)
  grid = _cellpoint_to_grid(points)
  pdat = _prepare_pointdata(pointdata)
  writevtk(grid,filebase,pointdata=pdat)
end

function _cellpoints_to_grid(points::CellPoints{D}) where D
  ps = Array{Point{D,Float64},1}(undef,(0,))
  p_to_cell = Array{Int,1}(undef,(0,))
  for (cell,p) in enumerate(points)
    for pj in p
      push!(ps,pj)
      push!(p_to_cell,cell)
    end
  end
  data, ptrs, ts, os = _prepare_cells(ps)
  grid = UnstructuredGrid(ps,data,ptrs,ts,os)
  (grid, p_to_cell)
end

function _cellpoint_to_grid(points::CellValue{Point{D,Float64}}) where D
  ps = collect(points)
  data, ptrs, ts, os = _prepare_cells(ps)
  UnstructuredGrid(ps,data,ptrs,ts,os)
end

function _prepare_cells(ps)
  data = [ i for i in 1:length(ps) ]
  ptrs = [ i for i in 1:(length(ps)+1) ]
  ts = ConstantCellValue( (), length(ps) )
  os = ConstantCellValue( 1, length(ps) )
  (data,ptrs,ts,os)
end

function _prepare_pointdata(pointdata)
  pdat = Dict()
  for (k,v) in pointdata
    pdat[k] = _prepare_data(v)
  end
  pdat
end

_prepare_data(v) = v

const IterData{T} = Union{CellValue{T},AbstractArray{T}}

function _prepare_data(v::IterData{<:VectorValue{D}}) where D
  a = collect(v)
  reshape(reinterpret(Float64,a),(D,length(a)))
end

function _prepare_data(v::IterData{<:VectorValue{2}})
  a = collect(v)
  b = reshape(reinterpret(Float64,a),(2,length(a)))
  z = zeros((1,size(b,2)))
  vcat(b,z)
end

function _prepare_data(v::IterData{<:TensorValue{D}}) where D
  a = collect(v)
  reshape(reinterpret(Float64,a),(D*D,length(a)))
end

_prepare_data(v::CellArray{<:Number}) = _flatten(v)

function _prepare_data(v::CellArray{<:VectorValue{D}}) where D
  a = _flatten(v)
  reshape(reinterpret(Float64,a),(D,length(a)))
end

function _prepare_data(v::CellArray{<:VectorValue{2}})
  a = _flatten(v)
  b = reshape(reinterpret(Float64,a),(2,length(a)))
  z = zeros((1,size(b,2)))
  vcat(b,z)
end

function _prepare_data(v::CellArray{<:TensorValue{D}}) where D
  a = _flatten(v)
  reshape(reinterpret(Float64,a),(D*D,length(a)))
end

struct VisualizationGrid{D,Z,G<:Grid{D,Z},C<:IndexCellValue{Int},P<:CellPoints{Z}} <: Grid{D,Z}
  grid::G
  coarsecells::C
  samplingpoints::P
end

points(vg::VisualizationGrid) = points(vg.grid)

cells(vg::VisualizationGrid) = cells(vg.grid)

celltypes(vg::VisualizationGrid) = celltypes(vg.grid)

function _writevtk(vg::VisualizationGrid,filebase,celldata,cellfields)
  cdata = _prepare_cdata(celldata,vg.coarsecells)
  pdata = _prepare_pdata(cellfields,vg.samplingpoints)
  _writevtk(vg.grid,filebase,cdata,pdata)
end

function _prepare_cdata(celldata,fine_to_coarse)
  cdata = Dict()
  for (k,v) in celldata
    acoarse = collect(v)
    afine = _allocate_afine(acoarse,length(fine_to_coarse))
    _fill_afine!(afine,acoarse,fine_to_coarse)
    cdata[k] = afine
  end
  k2 = "cellid"
  @assert ! haskey(cdata,k2)
  cdata[k2] = fine_to_coarse
  cdata
end

_allocate_afine(acoarse::Array{T},l) where T = Array{T,1}(undef,(l,))

function _fill_afine!(afine,acoarse,fine_to_coarse)
  for (i,coarse) in enumerate(fine_to_coarse)
    afine[i] = acoarse[coarse]
  end
end

function _prepare_pdata(cellfields,samplingpoints)
  pdata = Dict()
  for (k,v) in cellfields
    pdata[k] = _flatten(evaluate(v,samplingpoints))
  end
  pdata
end

function _writevtk(trian::Triangulation,filebase,nref,celldata,cellfields)
  vg = _visgrid(trian,nref)
  _writevtk(vg,filebase,celldata,cellfields)
end

function _visgrid(self::Triangulation,nref)
  grid, coarsecells, samplingpoints = _prepare_grid(CellRefFEs(self),CellGeomap(self),nref)
  VisualizationGrid(grid,coarsecells,samplingpoints)
end

function _refgrid end

function _refgrid(::Val{V},nref::Int) where V
  @notimplementedif nref != 0
  polytope = Polytope(V)
  Grid(polytope,length(V))
end

function _refgrid(::Val{(HEX_AXIS,)},nref::Int)
  n = 2^nref
  CartesianGrid(domain=(0.0,1.0),partition=(n,))
end

function _refgrid(::Val{(HEX_AXIS,HEX_AXIS)},nref::Int)
  n = 2^nref
  CartesianGrid(domain=(0.0,1.0,0.0,1.0),partition=(n,n))
end

function _refgrid(::Val{(HEX_AXIS,HEX_AXIS,HEX_AXIS)},nref::Int)
  n = 2^nref
  CartesianGrid(domain=(0.0,1.0,0.0,1.0,0.0,1.0),partition=(n,n,n))
end

function _prepare_grid(reffes::CellValue{<:RefFE},phi::CellGeomap{Z,D},nref) where {Z,D}
  @notimplemented
end

function _prepare_grid(reffes::ConstantCellValue{<:RefFE},phi::CellGeomap{Z,D},nref) where {Z,D}
  refgrid = _prepare_refgrid(reffes,nref)
  samplingpoints = _prepare_samplingpoints(refgrid,reffes)
  ps, offsets = _prepare_points(samplingpoints,points(refgrid),phi)
  data, ptrs, coarsecells = _prepare_cells(refgrid,offsets)
  ts = _prepare_celltypes(length(reffes),celltypes(refgrid))
  os = _prepare_cellorders(length(reffes),celltypes(refgrid))
  grid = UnstructuredGrid(ps,data,ptrs,ts,os)
  (grid, coarsecells, samplingpoints)
end

function _prepare_refgrid(reffes,nref)
  reffe = reffes.value
  poly = polytope(reffe)
  extrusion = poly.extrusion.array
  _refgrid(Val(Tuple(extrusion)),nref)
end

function _prepare_samplingpoints(refgrid,reffes)
  refpoints = collect(points(refgrid))
  ConstantCellValue(refpoints[:],length(reffes))
end

function _prepare_points(samplingpoints,refpoints,phi::CellGeomap{Z,D}) where {Z,D}
  xe = evaluate(phi,samplingpoints)
  offsets = Array{Int,1}(undef,(length(xe),))
  ps = Array{Point{D,Float64},1}(undef,(length(xe)*length(refpoints)))
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
  refextrusion = refcelltypes.value
  ConstantCellValue(refextrusion,ncells*length(refcelltypes) )
end

function _prepare_cellorders(ncells,refcelltypes::ConstantCellValue)
  ConstantCellValue(1,ncells*length(refcelltypes) )
end

function _flatten(ca::CellArray{T}) where T
  f = Vector{T}(undef,0)
  _flatten_kernel!(f,ca)
  f
end

function _flatten_kernel!(f,ca)
  for a in ca
    for ai in a
      push!(f,ai)
    end
  end
end

end # module Vtkio
