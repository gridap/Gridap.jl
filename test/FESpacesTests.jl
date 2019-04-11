using Numa, Test

using Numa.Quadratures
using Numa.Polytopes
using Numa.RefFEs
# using Numa.Meshes
using Numa.FESpaces

using Numa.Polytopes: PointInt
using Numa.CellValues

##
# Constant Pointer struct

struct ConstantPointerVector <: AbstractVector{Int}
  length::Int
  gap::Int
end

import Base: size
size(self::ConstantPointerVector) = (self.length,)

import Base:getindex
getindex(self::ConstantPointerVector, i::Int)::Int = (i-1)*self.gap+1
getindex(self::ConstantPointerVector, I::Vararg{Integer,1})::Int = (I[1]-1)*self.gap+1

import Base:setindex!
function setindex!(self::ConstantPointerVector, v, i)
    error("Index of ConstantPointerVector cannot be set")
end
##
point = ConstantPointerVector(10,3)
point[1]
point[3]
@test point == [(i-1)*point.gap+1 for i=1:length(point)]
##

##
data = [2,3,1,3,6,7,3,2,5,6,3,4]
ptrs = ConstantPointerVector(7,2)
celldata = CellVectorFromDataAndPtrs(data,ptrs)
@test celldata[end] == [3,4]
##

# Now return such structure from FESpace with same cell everywhere and conforming
# meshes


##
D = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
order=1
orders=order*ones(Int64,D)
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE(polytope,orders)
mesh = StructHexMesh(nparts)
fesp = FESpace(reffe,mesh)
@test fesp.l2giddof[end][end]==(nparts1d+1)^D
##

##
# Vector FE space
D = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
order=2
orders=order*ones(Int64,D)
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE(polytope,orders)
gps=[2,2]
quad=Quadratures.TensorProductQuadrature(orders = tuple(2*orders...))
shfscal = shfsps(reffe,quad.coords)
reffe = LagrangianRefFE(polytope,orders,2)
shftens = shfsps(reffe,quad.coords)
@test shftens[36,1,2,2]==shfscal[9,1,1]
##
points=quad.coords
reffe = LagrangianRefFE(polytope,orders,0)
gradscal = gradshfsps(reffe,points)
reffe = LagrangianRefFE(polytope,orders,1)
gradtens = gradshfsps(reffe,points)
@test gradtens[10,1,2,2]==gradscal[1,1,2]
@test gradtens[10,3,2,2]==gradscal[1,3,2]
@test gradtens[1,3,2,1]==gradscal[1,3,2]
@test gradtens[1,3,2,2]==0.0
@test gradtens[10,3,2,1]==0.0
##
