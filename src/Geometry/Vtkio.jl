
function writevtk(grid::Grid,filebase;celldata=Dict(),pointdata=Dict())
  points = vtkpoints(grid)
  cells = vtkcells(grid)
  vtkfile = vtk_grid(filebase, points, cells)
  for (k,v) in celldata
    vtk_cell_data(vtkfile, prepare_data(v), k)
  end
  for (k,v) in pointdata
    vtk_point_data(vtkfile, prepare_data(v), k)
  end
  outfiles = vtk_save(vtkfile)
end

function vtkpoints(grid::Grid{D}) where D
  x = points(grid)
  xflat = collect(x)
  reshape(reinterpret(Float64,xflat),(D,length(x)))
end

# @fverdugo this allocates a lot of small objects
# Not very crucial since it is for visualization
# but it would be nice to have a better way
function vtkcells(grid::Grid)
  types = vtkcelltypedict()
  nodes = vtkcellnodesdict()
  c = celltypes(grid)
  n = cells(grid)
  [ MeshCell(types[encode_extrusion(ci)], ni[nodes[encode_extrusion(ci)]])
     for (ci,ni) in zip(c,n) ] 
end

"""
Generates the lookup table (as a Dict) in order to convert between
Numa Polytope identifiers into VTK cell type identifiers
"""
function vtkcelltypedict()
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

"""
Generates the lookup table (as a Dict) in order to convert between
Numa Polytope corner numbering into VTK corner numbering
"""
function vtkcellnodesdict()
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

function writevtk(points::CellPoints,filebase;celldata=Dict(),pointdata=Dict())
  grid, p_to_cell = cellpoints_to_grid(points)
  pdat = prepare_pointdata(pointdata)
  k = "cellid"
  @assert ! haskey(pdat,k)
  pdat[k] = p_to_cell
  writevtk(grid,filebase,pointdata=pdat)
end

function writevtk(points::CellValue{Point{D}} where D,filebase;celldata=Dict(),pointdata=Dict())
  grid = cellpoint_to_grid(points)
  pdat = prepare_pointdata(pointdata)
  writevtk(grid,filebase,pointdata=pdat)
end

function cellpoints_to_grid(points::CellPoints{D}) where D
  ps = Array{Point{D},1}(undef,(0,))
  p_to_cell = Array{Int,1}(undef,(0,))
  for (cell,p) in enumerate(points)
    for pj in p
      push!(ps,pj)
      push!(p_to_cell,cell)
    end
  end
  data, ptrs, ts = prepare_cells(ps)
  grid = UnstructuredGrid(ps,data,ptrs,ts)
  (grid, p_to_cell)
end

function cellpoint_to_grid(points::CellValue{Point{D}}) where D
  ps = collect(points)
  data, ptrs, ts = prepare_cells(ps)
  UnstructuredGrid(ps,data,ptrs,ts)
end

function prepare_cells(ps)
  data = [ i for i in 1:length(ps) ]
  ptrs = [ i for i in 1:(length(ps)+1) ]
  ts = [ () for i in 1:length(ps) ]
  (data,ptrs,ts)
end

function prepare_pointdata(pointdata)
  pdat = Dict()
  for (k,v) in pointdata
    pdat[k] = prepare_data(v)
  end
  pdat
end

prepare_data(v) = v

function prepare_data(v::IterData{<:VectorValue{D}}) where D
  a = collect(v)
  reshape(reinterpret(Float64,a),(D,length(a)))
end

function prepare_data(v::IterData{<:VectorValue{2}})
  a = collect(v)
  b = reshape(reinterpret(Float64,a),(2,length(a)))
  z = zeros((1,size(b,2)))
  vcat(b,z)
end

function prepare_data(v::IterData{<:TensorValue{D}}) where D
  a = collect(v)
  reshape(reinterpret(Float64,a),(D*D,length(a)))
end

prepare_data(v::CellArray{<:Number}) = collect(flatten(v))

function prepare_data(v::CellArray{<:VectorValue{D}}) where D
  a = collect(flatten(v))
  reshape(reinterpret(Float64,a),(D,length(a)))
end

function prepare_data(v::CellArray{<:VectorValue{2}})
  a = collect(flatten(v))
  b = reshape(reinterpret(Float64,a),(2,length(a)))
  z = zeros((1,size(b,2)))
  vcat(b,z)
end

function prepare_data(v::CellArray{<:TensorValue{D}}) where D
  a = collect(flatten(v))
  reshape(reinterpret(Float64,a),(D*D,length(a)))
end

