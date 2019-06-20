module Triangulations

using Test
using Gridap
using Gridap.Helpers
using StaticArrays

export Triangulation
export CellRefFEs
export CartesianTriangulation
export test_triangulation
export ncells
import Gridap: CellQuadrature
import Gridap: CellPoints
import Gridap: CellBasis
import Gridap: CellGeomap
import Gridap: CellField
import Base: IndexStyle
import Base: size
import Base: getindex

"""
Minimal interface for a mesh used for numerical integration
I.e., cellwise nodal coordinates and notion about interpolation of these
cell coordinates within the cells.
"""
abstract type Triangulation{Z,D} end

function CellPoints(::Triangulation{Z,D})::CellPoints{D,Float64} where {Z,D}
 @abstractmethod
end

function CellRefFEs(
  trian::Triangulation{Z})::CellValue{LagrangianRefFE{Z,Float64}} where Z
  @abstractmethod
end

function CellBasis(trian::Triangulation{Z})::CellBasis{Z,Float64} where Z
 @abstractmethod
end

function CellGeomap(self::Triangulation)
  coords = CellPoints(self)
  basis = CellBasis(self)
  lincomb(basis,coords)
end

function ncells(self::Triangulation)
  coords = CellPoints(self)
  length(coords)
end

# Testers

function test_triangulation(trian::Triangulation{Z,D}) where {Z,D}
  basis = CellBasis(trian)
  @test isa(basis,CellBasis{Z,Float64})
  coords = CellPoints(trian)
  @test isa(coords,CellPoints{D,Float64})
  @test ncells(trian) == length(coords)
  @test ncells(trian) == length(basis)
  phi = CellGeomap(trian)
  @test isa(phi,CellGeomap{Z,D,Float64})
  jac = gradient(phi)
  @test isa(jac,CellField{D,TensorValue{D,Float64,D*D}})
  reffes = CellRefFEs(trian)
  @test isa(reffes,CellValue{LagrangianRefFE{Z,Float64}})
end

# Factories

function CellField(trian::Triangulation,fun::Function)
  phi = CellGeomap(trian)
  compose(fun,phi)
end

function CellField(
  trian::Triangulation{D,Z}, fun::Function, u::Vararg{<:CellField{Z}}) where {D,Z}
  phi = CellGeomap(trian)
  compose(fun,phi,u...)
end

function CellBasis(
  trian::Triangulation{D,Z},
  fun::Function,
  b::CellBasis{Z},
  u::Vararg{<:CellField{Z}}) where {D,Z}
  phi = CellGeomap(trian)
  compose(fun,phi,b,u...)
end

"""
Factory function to create CellQuadrature objects in a convenient way
"""
function CellQuadrature(trian::Triangulation;order::Int)
  _quadrature(CellRefFEs(trian),order)
end

_quadrature(reffes,order) = @notimplemented

function _quadrature(reffes::ConstantCellValue,order)
  reffe = reffes.value
  poly = polytope(reffe)
  t = poly.extrusion.data
  q = Quadrature(t,order=order)
  ncells = length(reffes)
  ConstantCellValue(q,ncells)
  ConstantCellQuadrature(q,ncells)
end

# Concrete implementations

struct CartesianTriangulation{D,C,B,R} <: Triangulation{D,D}
  coords::C
  basis::B
  reffes::R
end

function CartesianTriangulation(partition,domain=nothing)
  _p = _setup_partition(partition)
  _d = _setup_domain(length(_p),domain)
  D = length(_p)
  coords = CartesianCellPoints(_p,_d)
  ncells = prod(_p)
  basis, reffes = _setup_basis(D,ncells)
  C = typeof(coords)
  B = typeof(basis)
  R = typeof(reffes)
  CartesianTriangulation{D,C,B,R}(coords,basis,reffes)
end

CellBasis(trian::CartesianTriangulation) = trian.basis

CellPoints(trian::CartesianTriangulation) = trian.coords

CellRefFEs(trian::CartesianTriangulation) = trian.reffes

_setup_partition(partition::Vector{Int}) = tuple(partition...)

_setup_partition(partition::NTuple{D,Int} where D) = partition

function _setup_domain(D,::Nothing)
  d =  [ i*(-1)^j for i in ones(D) for j in 1:2 ]
  _setup_domain(D,d)
end

_setup_domain(D,domain::Vector{Float64}) = tuple(domain...)

_setup_domain(D,domain::NTuple{D,Int} where D) = domain

function _setup_basis(D,ncells)
  code = @SVector fill(HEX_AXIS,D)
  polytope = Polytope(code)
  order = 1
  orders = fill(order,D)
  reffe = LagrangianRefFE{D,Float64}(polytope,orders)
  basis = shfbasis(reffe)
  b = ConstantCellMap(basis, ncells)
  r = ConstantCellValue(reffe,ncells)
  (b,r)
end

struct CartesianCellPoints{D} <: IndexCellArray{Point{D,Float64},1,Vector{Point{D,Float64}},D}
  p0::Point{D,Float64}
  dp::Point{D,Float64}
  partition::NTuple{D,Int}
  v::Vector{Point{D,Float64}}
end

function CartesianCellPoints(
  partition::NTuple{D,Int}, domain::NTuple{N,Float64}) where {D,N}

  dim_to_limits = tuple([(domain[2*i-1],domain[2*i]) for i in 1:D ]...)
  p0 = zero(MVector{D,Float64})
  dp = zero(MVector{D,Float64})
  for d in 1:D
    p0[d] = dim_to_limits[d][1]
    dp[d] = (dim_to_limits[d][2] - dim_to_limits[d][1]) / partition[d]
  end
  v = zeros(Point{D,Float64},2^D)
  CartesianCellPoints(Point(p0),Point(dp),partition,v)
end

IndexStyle(::Type{CartesianCellPoints{D}}) where D = IndexCartesian()

size(x::CartesianCellPoints) = x.partition

function getindex(self::CartesianCellPoints{D}, cell::Vararg{<:Integer,D}) where D
  lsize = @SVector fill(2,D)
  cis = CartesianIndices(lsize.data)
  p = zero(MVector{D,Float64})
  p0 = self.p0
  dp = self.dp
  for (k,ci) in enumerate(cis)
    for d in 1:D
      p[d] = p0[d] + ( (cell[d]-1) + (ci[d]-1) )*dp[d]
    end
    self.v[k] = p
  end
  self.v
end

end # module
