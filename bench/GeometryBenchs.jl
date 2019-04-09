module GeometryBenchs

using Numa.Geometry: CartesianGridCoords
using Numa.Geometry: CartesianGridConnectivity

function doloop(a)
  for ai in a
  end
end

l = 1000000

println("+++ GeometryBenchs ( length = $l ) +++")

domain = ((0.0,1.0),(-1.0,2.0),(2.0,3.0))
npoints = (101,101,101)
x = CartesianGridCoords(domain,npoints)

print("CartesianGridCoords ->"); @time doloop(x)
print("CartesianGridCoords ->"); @time doloop(x)

ncells = (100,100,100)
t = CartesianGridConnectivity(ncells)

print("CartesianGridConnectivity ->"); @time doloop(t)
print("CartesianGridConnectivity ->"); @time doloop(t)

end # module GeometryBenchs
