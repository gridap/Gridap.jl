##
using Numa, Test

using Numa.Quadratures
using Numa.Polytopes
using Numa.RefFEs
# using Numa.Meshes
using Numa.FESpaces
using Numa.FESpaces: ConformingFESpace

using Numa.Polytopes: PointInt
using Numa.CellValues

using Numa.Maps
using Numa.FieldValues

using Numa.Meshes

using Numa.CellValues: CellVectorByComposition

using Numa.CellMaps

using Numa.Geometry
using Numa.Geometry.Cartesian

import Numa: gradient, ∇

using Numa.CellIntegration

D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = triangulation(grid) # Generates the Triangulation associated with this grid
graph = gridgraph(grid) # Generates the GridGraph associated with this grid.
phi = geomap(trian)

order=1
orders=order*ones(Int64,D)
pol_array = celltypes(trian)
extrusion = PointInt{D}(pol_array[1])
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
basis = reffe.shfbasis
cellb = ConstantCellValue(basis, ncells(trian))
quad = quadrature(trian,order=2)
gps = coordinates(quad)

fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
# fun(x::Point{2}) = x[1]^2
# gradfun(x::Point{2}) = VectorValue(2*x[1], 0.0)
# gradient(::typeof(fun)) = gradfun

fesp = ConformingFESpace(reffe,trian,graph)
assembler = ConformingAssembler(fesp)
funh = interpolate(fun, fesp)

p_arr = []
for pi in points(grid)
  global p_arr
  p_arr = [p_arr..., pi]
end
p_arr
fun.(p_arr)

@test funh.coeffs.gid_to_val == fun.(p_arr)

# Next, put negative values, etc...
D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
trian = triangulation(grid) # Generates the Triangulation associated with this grid
quad = quadrature(trian,order=4)
fun(x::Point{2}) = 1.0
ksca = integrate(fun, trian, quad)
int = sum(ksca)
@test int ≈ 1.0
# @santiagobadia : I would say that this is not the result we want...
# we want the whole integral over the whole domain ... 1.0

physbasis = attachgeomap(cellb,phi);
ab(v,u) = inner(∇(v),∇(u)) #+ inner(v,u)
am(v,u) = inner(v,u)
V = physbasis;
U = physbasis;

uphys = fun ∘ phi

evaluate(phi, gps)
evaluate(uphys, gps)
ksca = integrate(ab(uphys,uphys),trian,quad)
int = sum(ksca)
@test int ≈ 1.0
##

kvec = integrate(ab(V,uphys),trian,quad)
# @santiagobadia : Problem when showing this result, in iterate
kvec = integrate(am(V,uphys),trian,quad)
typeof(kvec)
kmat = integrate(ab(V,U),trian,quad)
# @santiagobadia : Problem in evaluation...

# using Numa.CellMaps: CellFieldFromCompose
# function gradient(self::CellFieldFromCompose)
#   gradop = gradient(self.op)
#   CellFieldFromCompose(gradop,self.a)
# end


gradient(fun)
