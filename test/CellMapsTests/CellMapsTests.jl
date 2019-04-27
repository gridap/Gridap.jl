##
using Numa
using Test

using Numa.Quadratures
using Numa.CellQuadratures
import Numa.CellIntegration: cellcoordinates, cellbasis

using Numa.CellValues
using Numa.CellValues: IndexCellArray
import Numa.CellValues: cellsize

using Numa.Polytopes
using Numa.Polytopes: PointInt

using Numa.RefFEs
using Numa.FieldValues

import Numa: gradient
using Numa.CellValues: ConstantCellValue

using Numa.FESpaces: ConformingFESpace

using Numa.Maps
using Numa.Maps: AnalyticalField

using Numa.Geometry
using Numa.Geometry.Cartesian
using Numa.Geometry.Unstructured

using Numa.CellMaps

D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
grid = CartesianGrid(partition=nparts_t,domain=(0,1,0,1),order=1) # domain, and order are optional
grid = UnstructuredGrid(grid)
trian = triangulation(grid) # Generates the Triangulation associated with this grid
# graph = gridgraph(grid)
## # Generates the GridGraph associated with this grid.
phi = geomap(trian)
using Numa.CellMaps: CellFieldFromExpand
typeof(phi) <: CellFieldFromExpand
l = prod(nparts_t)
refquad = TensorProductQuadrature(orders=(5,4))
refpoints = coordinates(refquad)
quad = ConstantCellQuadrature(refquad,l)
quad
p = coordinates(quad)
x = evaluate(phi,p)
##
fun(x::Point{2}) = x[2]
gradfun(x::Point{2}) = VectorValue(0.0, 1.0)
gradient(::typeof(fun)) = gradfun
f = AnalyticalField(fun,2)
using Numa.Maps
using Numa.CellMaps: ConstantCellMap
@test typeof(f) <: Map{Point{D},1,Float64,1}
ccm = ConstantCellMap(f,10)
@test typeof(ccm[1]) <: AnalyticalField
@test length(ccm) ==  10
@test size(ccm) == (10,)
res = evaluate(ccm,x)
##
AnalyticalField{2,Float64,typeof(fun)} <: Map{Point{D},1,Float64,1}
isa(ccm,ConstantCellValue{AnalyticalField{2,Float64,typeof(fun)}})
isa(ccm,ConstantCellMap{Point{D},1,Float64,1})














































#
# l = 10
#
# refquad = TensorProductQuadrature(orders=(5,4))
# refpoints = coordinates(refquad)
#
# quad = ConstantCellQuadrature(refquad,l)
# points = coordinates(quad)
#
# polytope = Polytope(Polytopes.PointInt{2}(1,1))
# reffe = LagrangianRefFE{2,ScalarValue}(polytope,[1,1])
# refbasis = reffe.shfbasis
#
# refvals = evaluate(refbasis,refpoints)
# ##
# import Numa.CellMaps: IterConstantCellMapValues
# vals = IterConstantCellMapValues(refbasis,points)
# typeof(vals)
# @test isa(vals,CellBasisValues)
#
# for refvals2 in vals
#   @assert refvals2 == refvals
# end
#
# basis = ConstantCellMap(refbasis, length(points))
#
# vals = evaluate(basis,points)
#
#
# for refvals2 in vals
#   @assert refvals2 == refvals
# end
#
# refbasisgrad = gradient(refbasis)
#
# refvalsgrad = evaluate(refbasisgrad,refpoints)
#
# basisgrad = gradient(basis)
#
# valsgrad = evaluate(basisgrad,points)
#
# @test isa(valsgrad,CellBasisValues{VectorValue{2}})
#
#
# for refvalsgrad2 in valsgrad
#   @assert refvalsgrad2 == refvalsgrad
# end
#
# valsgrad = IterConstantCellMapValues(refbasisgrad,points)
#
# @test isa(valsgrad,CellBasisValues{VectorValue{2}})
#
# for refvalsgrad2 in valsgrad
#   @assert refvalsgrad2 == refvalsgrad
# end
# ##
#
#
# @eval begin
#
#   ufun(x::Point{2}) = x[1]*x[2] + x[1]
#
#   gradufun(x::Point{2}) = VectorValue(x[2]+1.0,x[1])
#
#   gradient(::typeof(ufun)) = gradufun
#
# end
#
# @test isa(phi,CellField)
# @test isa(ufun,Function)
#
# cfield = compose(ufun,phi)
#
# @test isa(cfield,CellField{2,Float64})
#
# cfieldgrad = gradient(cfield)
#
# @test isa(cfieldgrad,CellField{2,VectorValue{2}})
#
# uatx = evaluate(cfield,points)
#
# ugradatx = evaluate(cfieldgrad,points)
#
# x = evaluate(phi,points)
#
# for (ui,uigrad,xi) in zip(uatx,ugradatx,x)
#   @assert ui == ufun.(xi)
#   @assert uigrad == gradufun.(xi)
# end
