##
using Numa
using Test

using Numa.Quadratures
using Numa.CellQuadratures
import Numa.CellIntegration: cellcoordinates, cellbasis

using Numa.CellValues: IndexCellArray
import Numa.CellValues: cellsize

using Numa.Polytopes
using Numa.Polytopes: PointInt

using Numa.RefFEs
using Numa.FieldValues

import Numa: gradient
using Numa.CellValues: ConstantCellValue

using Numa.FESpaces: ConformingFESpace

using Numa.Maps: AnalyticalField

using Numa.Geometry
using Numa.Geometry.Cartesian

using Numa.CellMaps

D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = triangulation(grid) # Generates the Triangulation associated with this grid
graph = gridgraph(grid) # Generates the GridGraph associated with this grid.
phi = geomap(trian)
fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
f = AnalyticalField(fun,2)

using Numa.Maps
using Numa.CellMaps: ConstantCellMap
@test typeof(f) <: Map{Point{D},1,Float64,1}
ccm = ConstantCellMap(f,10)
@test typeof(ccm[1]) <: AnalyticalField
@test length(ccm) ==  10
@test size(ccm) == (10,)

p = Point{D}(2,2)
@test ccm[1].fun(p) == 2.0

quad = quadrature(trian,order=2)
points = quad.coords.array
@test ccm[1].fun.(points)[1]  ≈ - 1*sqrt(3)/3
using Numa.CellValues
points = quad.coords
@test typeof(points) <: CellArray{Point{D},1}
@test typeof(f) <: Map{Point{D},1,Float64,1}
typeof(ccm)
val = evaluate(ccm,points)
@test cellsize(points) == cellsize(val)
##
