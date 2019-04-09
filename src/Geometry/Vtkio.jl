
function writevtk(grid::Grid,filebase)
  points = vtkpoints(grid)
  cells = vtkcells(grid)
  vtkfile = vtk_grid(filebase, points, cells)
  outfiles = vtk_save(vtkfile)
  reportfiles(outfiles)
end

function vtkpoints(grid::Grid{D}) where D
  x = coordinates(grid)
  points = collect(x)
  reshape(reinterpret(Float64,points),(D,length(x)))
end

function vtkcells(grid::Grid)
  types = vtkcelltypedict()
  nodes = vtkcellnodesdict()
  c = celltypes(grid)
  n = connectivity(grid)
  [  MeshCell(types[ci], ni[nodes[ci]]) for (ci,ni) in zip(c,n) ] 
end

function reportfiles(outfiles)
  for file in outfiles
    println("Generated results file: $file")
  end
end

function vtkcelltypedict()
  d = Dict()
  h = HEX_AXIS
  t = TET_AXIS
  d[()] = VTK_VERTEX
  d[(t,)] = VTK_LINE
  d[(h,)] = VTK_LINE
  d[(t,t)] = VTK_TRIANGLE
  d[(h,h)] = VTK_QUAD
  d[(t,t,t)] = VTK_TETRA
  d[(h,h,h)] = VTK_HEXAHEDRON
  d
end

function vtkcellnodesdict()
  d = Dict()
  h = HEX_AXIS
  t = TET_AXIS
  d[()] = [1,]
  d[(t,)] = [1,2]
  d[(h,)] = d[(t,)]
  d[(t,t)] = [1,2,3]
  d[(h,h)] = [1,2,4,3]
  d[(t,t,t)] = [1,2,3,4]
  d[(h,h,h)] = [1,2,4,3,5,6,8,7]
  d
end
