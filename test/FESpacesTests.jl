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

using Numa.Fields
using Numa.FieldValues

using Numa.Meshes

##
D = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
order=3
orders=order*ones(Int64,D)
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
mesh = StructHexMesh(nparts)
fesp = ConformingFESpace{D,D,ScalarValue}(reffe,mesh)
using Numa.FESpaces: globaldofs
gldofs = globaldofs(fesp)
@test gldofs[end][end]==(nparts1d*order+1)^D
##


gldofs
mesh.cellvefs

# @santiagobadia : I want to compose two cellvectors... to get all cellvefs
# put the result in scratch data...

using Numa.CellValues: IndexCellArray
using Numa.CellValues: IndexCellVector
using Numa.CellValues: CachedArray

typeof(gldofs) <: IndexCellVector{Int64}
typeof(mesh.cellvefs) <: IndexCellVector{Int64}

struct CellVectorByComposition{T,L<:IndexCellVector{Int},V<:IndexCellVector{T}} <: IndexCellArray{T,1,CachedArray{T,1,Array{T,1}},1}
  cell_to_x::L
  x_to_vals::V
  cv::CachedArray{T,1,Array{T,1}}
end

function CellVectorByComposition(cell_to_x::IndexCellVector{Int}, x_to_vals::IndexCellVector{T}) where T
  L = typeof(cell_to_x)
  V = typeof(x_to_vals)
  # Here something far more complicated
  a = Vector{T}(undef,(celllength(cell_to_x),))
  cv = CachedArray(a,size(a))
  CellVectorByComposition{T,L,V}(cell_to_x, x_to_vals, cv)
end

celldofs = CellVectorByComposition(mesh.cellvefs,gldofs)

using Base: @propagate_inbounds

@propagate_inbounds function getindex(self::CellVectorFromLocalToGlobal,cell::Int)
  cell_to_x = self.cell_to_x[cell]
  l = 0
  for x in cell_to_x
    l += length(self.x_to_vals[x])
  end
  setsize!(self.cv,(l,))
  l = 1
  for x in cell_to_x
    for val in self.x_to_vals[x]
      self.cv[l] = val
      l += 1
    end
  end
  self.cv
end

import Base: size, IndexStyle
size(self::CellVectorFromLocalToGlobal) = (length(self.cell_to_x),)

IndexStyle(::Type{CellVectorFromLocalToGlobal{T,L,V}}) where {T,L,V} = IndexLinear()

cellsize(self::CellVectorFromLocalToGlobal) = cellsize(self.cell_to_x)*cellsize(self.x_to_vals)





























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
