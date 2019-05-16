module GeometryBenchs

using Gridap.Geometry
using Gridap.Geometry.Cartesian
using Gridap.Geometry.Unstructured

function doloop(a)
  for ai in a
  end
end

l = 1000000

println("+++ GeometryBenchs ( length = $l ) +++")

grid = CartesianGrid(domain=(0.0,1.0,-1.0,2.0,2.0,3.0),partition=(100,100,100))

x = points(grid)

print("CartesianGridPoints ->"); @time doloop(x)
print("CartesianGridPoints ->"); @time doloop(x)

t = cells(grid)

print("CartesianGridCells ->"); @time doloop(t)
print("CartesianGridCells ->"); @time doloop(t)

print("FlexibleUnstructuredConstructor ->"); @time ugrid = FlexibleUnstructuredGrid(grid)
print("FlexibleUnstructuredConstructor ->"); @time ugrid = FlexibleUnstructuredGrid(grid)

x = points(ugrid)

print("FlexibleUnstructuredGridPoints ->"); @time doloop(x)
print("FlexibleUnstructuredGridPoints ->"); @time doloop(x)

t = cells(ugrid)

print("FlexibleUnstructuredGridCells ->"); @time doloop(t)
print("FlexibleUnstructuredGridCells ->"); @time doloop(t)

print("UnstructuredConstructor ->"); @time ugrid = UnstructuredGrid(grid)
print("UnstructuredConstructor ->"); @time ugrid = UnstructuredGrid(grid)

x = points(ugrid)

print("UnstructuredGridPoints ->"); @time doloop(x)
print("UnstructuredGridPoints ->"); @time doloop(x)

t = cells(ugrid)

print("UnstructuredGridCells ->"); @time doloop(t)
print("UnstructuredGridCells ->"); @time doloop(t)

end # module GeometryBenchs
