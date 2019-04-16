##
using Numa
using Numa.Quadratures
using Numa.CellQuadratures
# using Numa.CellIntegration
# using Numa.CellValues
using Numa.CellFunctions
import Numa.CellIntegration: cellcoordinates, cellbasis

using Numa.CellValues: IndexCellArray
import Numa.CellValues: cellsize
using Numa.Polytopes
using Numa.Polytopes: PointInt
using Numa.RefFEs
using Numa.FieldValues

import Numa: gradient
using Numa.CellValues: ConstantCellValue
include("CellIntegrationTestsMocks.jl")

using Numa.Meshes
using Numa.FESpaces: ConformingFESpace

##
##
D=2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
nparts_t = tuple(nparts...)
order=1
orders=order*ones(Int64,D)
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytopes.Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope,orders)
basis = reffe.shfbasis
cellb = CellBasisFromSingleInterpolation(basis)
imesh = DummyIntegrationMesh2D(partition=nparts_t)
refquad = TensorProductQuadrature(orders=(2,2))
meshcoords = cellcoordinates(imesh)
ncells = length(meshcoords)
quad = ConstantCellQuadrature(refquad,ncells)
phi = geomap(imesh)
basis = cellbasis(imesh)
physbasis = attachgeomap(basis,phi)
ab(v,u) = inner(∇(v),∇(u)) #+ inner(v,u)
V = physbasis
U = physbasis
# fun(x::Point{2}) = x[1]*x[2] + x[1]
# gradfun(x::Point{2}) = VectorValue(x[2] + 1.0, x[1])
fun(x::Point{2}) = x[1]
gradfun(x::Point{2}) = VectorValue(1.0, 0.0)
gradient(::typeof(fun)) = gradfun
uphys = fun ∘ phi


















##


"""
Concrete implementation of CellFunction for the case of the same function on all cells
"""
struct ConstantCellFunction{S,M,T,N}
  A::Function{S,M,T,N}  # Or put an acceptable name here
end

function CellBasisFromSingleInterpolation(f::ConstantCellFunction{S,M,T,N}) where {S,M,T,N}
  ConstantCellFunction{S,M,T,N}(f)
end


"""
Cell-wise field created from a `Field`
"""
struct CellFieldFromField{D,T,F<:Field{D,T}} <: CellField{D,T}
  field::F
end

function evaluate(self::CellFieldFromField{D,T},points::CellPoints{D}) where D
  CellFieldValuesFromField(self.field,points)
end

function gradient(self::CellFieldFromField)
  gradfield = gradient(self.field)
  CellFieldFromField(gradfield)
end

"""
Result of the lazy evaluation of a `CellFieldFromField` on a cell-wise array
of points `CellPoints{D}`
"""
struct CellFieldValuesFromField{D,T,F,C} # <: IndexCellArray{T,1,C}
  # @santiagobadia : Not sure from where to extend
  field::F
  cellpoints::C
  cv::CachedVector{T,Vector{T}}
end

function CellFieldValuesFromField(
  field::Field{D,T},cellpoints::CellPoints{D}) where {D,T}
  F = typeof(field)
  C = typeof(cellpoints)
  a = Vector{T}(undef,celllength(cellpoints))
  cv = CachedArray(a)
  CellFieldValuesFromField{D,T,F,C}(field,cellpoints,cv)
end

# @santiagobadia : I am probably missing interface methods

@propagate_inbounds function getindex(self::CellFieldValuesFromField,cell::Int)
  points = self.cellpoints[cell]
  setsize!(self.cv,(length(points),))
  evaluate!(self.field,points,cv)
end

inputcellarray(self::CellFieldValuesFromField) = self.cellpoints

computesize(self::CellFieldValuesFromField, asize) = asize[1]

function computevals!(self::CellFieldValuesFromField, a, v)
  evaluate!(self.field,a,v)
end
