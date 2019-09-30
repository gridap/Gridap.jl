module FESpaces

using Test
using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery

export free_dofs
export diri_dofs
export FESpace
export test_fe_space
export TestFESpace
export TrialFESpace
export num_free_dofs
export num_diri_dofs
export diri_tags
export apply_constraints
export apply_constraints_rows
export apply_constraints_cols
export celldofids
export interpolate_values
export interpolate_diri_values
export value_type
export FESpaceWithDirichletData
import Gridap: CellField
import Gridap: CellBasis
import Gridap: Triangulation

# Interfaces

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank and entry type).
A FE space has to be understood as the direct sum of two spaces, the one with
arbitrary free values and zero Dirichlet data and the one with zero free values
and a given Dirichlet data
"""
abstract type FESpace{D,Z,T} end

num_free_dofs(::FESpace)::Int = @abstractmethod

num_diri_dofs(::FESpace)::Int = @abstractmethod

diri_tags(::FESpace)::Vector{Int} = @abstractmethod

function apply_constraints(::FESpace, cellvec::CellVector, cellids::CellNumber)::CellVector
  @abstractmethod
end

function apply_constraints_rows(::FESpace, cellmat::CellMatrix, cellids::CellNumber)::CellMatrix
  @abstractmethod
end

function apply_constraints_cols(::FESpace, cellmat::CellMatrix, cellids::CellNumber)::CellMatrix
  @abstractmethod
end

"""
Cell DOFs ids after applying constraints
"""
function celldofids(::FESpace)::CellVector{Int}
  @abstractmethod
end

"""
Interpolates a function and returns a vector of free values and another vector
of dirichlet ones. It is templatized with respect to E, the type of field
value in the `FESpace{D,Z,T}` at hand, i.e.,  `E = eltype(T)`
"""
function interpolate_values(::FESpace,::Function)::Tuple{Vector{E},Vector{E}} where E
  @abstractmethod
end

function interpolate_values(::FESpace{D,Z,T},::T) where {D,Z,T}
  @abstractmethod
end

"""
Given a FESpace and its array of labels that represent the Dirichlet boundary
in the geometrical model and an array of functions for every Dirichlet
boundary label, this method interpolates all these functions on their
boundaries, and provides the global vector of Dirichlet DOFs
"""
function interpolate_diri_values(::FESpace, funs::Vector{<:Function})::Vector{E} where E
  @abstractmethod
end

function interpolate_diri_values( ::FESpace{D,Z,T}, ::Vector{T}) where {D,Z,T}
  @abstractmethod
end

function interpolate_diri_values(this::FESpace, fun::Function)
  tags = diri_tags(this)
  interpolate_diri_values(this,fill(fun,length(tags)))
end

function interpolate_diri_values(this::FESpace{D,Z,T}, val::T) where {D,Z,T}
  tags = diri_tags(this)
  interpolate_diri_values(this,fill(val,length(tags)))
end

"""
Returns the CellField that represents the FE function thorough its free and
dirichlet values. E = eltype(T)
"""
function CellField(
  ::FESpace{D,Z,T},free_dofs::AbstractVector{E},diri_dofs::AbstractVector{E})::CellField{Z,T} where {D,Z,T,E}
  @abstractmethod
end

function CellBasis(::FESpace{D,Z,T})::CellBasis{Z,T} where {D,Z,T}
  @abstractmethod
end

function Triangulation(::FESpace{D,Z})::Triangulation{Z,D} where {Z,D}
  @abstractmethod
end

value_type(::FESpace{D,Z,T}) where {D,Z,T} = T

function TestFESpace(this::FESpace{D,Z,T}) where {D,Z,T}
  E = eltype(T)
  dv = zeros(E,num_diri_dofs(this))
  return FESpaceWithDirichletData(this, dv)
end

function TrialFESpace( this::FESpace, funs::Vector{<:Function}=Function[])
  _TrialFESpace(this,funs)
end

function TrialFESpace( this::FESpace{D,Z,T}, vals::Vector{T}) where {D,Z,T}
  _TrialFESpace(this,vals)
end

function _TrialFESpace( this, funs)
  tags = diri_tags(this)
  @assert length(tags) == length(funs)
  dv = interpolate_diri_values(this,funs)
  return FESpaceWithDirichletData(this, dv)
end

function TrialFESpace( this::FESpace, fun::Function)
  dv = interpolate_diri_values(this,fun)
  return FESpaceWithDirichletData(this, dv)
end

function TrialFESpace( this::FESpace{D,Z,T}, val::T) where {D,Z,T}
  dv = interpolate_diri_values(this,val)
  return FESpaceWithDirichletData(this, dv)
end

# Testers

function test_fe_space(
  fespace::FESpace{D,Z,T},
  nfree::Integer,
  ndiri::Integer,
  cellmat::CellMatrix,
  cellvec::CellVector,
  ufun::Function) where {D,Z,T}

  @test num_free_dofs(fespace) == nfree
  @test num_diri_dofs(fespace) == ndiri
  tags = diri_tags(fespace)
  @test isa(tags,Vector{Int})

  basis = CellBasis(fespace)
  @test isa(basis,CellBasis)

  nc = length(basis)
  cellids = IdentityCellNumber(Int,nc)

  cv = apply_constraints(fespace,cellvec,cellids)
  @test isa(cv,CellVector)
  cm = apply_constraints_rows(fespace,cellmat,cellids)
  @test isa(cm,CellMatrix)
  cm = apply_constraints_cols(fespace,cellmat,cellids)
  @test isa(cm,CellMatrix)

  cell_to_dofs = celldofids(fespace)
  @test isa(cell_to_dofs,CellVector{<:Integer})

  freevals, dirivals = interpolate_values(fespace,ufun)
  @test isa(freevals,AbstractVector)
  @test isa(dirivals,AbstractVector)
  @test length(freevals) == nfree
  @test length(dirivals) == ndiri

  z = zero(T)
  freevals, dirivals = interpolate_values(fespace,z)
  @test isa(freevals,AbstractVector)
  @test isa(dirivals,AbstractVector)
  @test length(freevals) == nfree
  @test length(dirivals) == ndiri

  dirivals = interpolate_diri_values(fespace,ufun)
  @test isa(dirivals,AbstractVector)
  @test length(dirivals) == ndiri

  cf = CellField(fespace,freevals,dirivals)
  @test isa(cf,CellField)

  trian = Triangulation(fespace)
  @test isa(trian,Triangulation{Z,D})

  dirivals = interpolate_diri_values(fespace,z)
  @test isa(dirivals,AbstractVector)
  @test length(dirivals) == ndiri

end

# Pretty printing

import Base: show

function show(io::IO,self::FESpace{D,Z,T}) where {D,Z,T}
  print(io,"$(nameof(typeof(self))) object")
end

function show(io::IO,::MIME"text/plain",self::FESpace{D,Z,T}) where {D,Z,T}
  show(io,self)
  print(io,":")
  print(io,"\n physdim: $D")
  print(io,"\n refdim: $Z")
  print(io,"\n valuetype: $T")
end

# Helpers

"""
FESpace whose Dirichlet component has been constrained
"""
struct FESpaceWithDirichletData{D,Z,T,E,V<:FESpace{D,Z,T}} <: FESpace{D,Z,T}
  fespace::V
  diri_dofs::Vector{E}
end

function free_dofs end

diri_dofs(f::FESpaceWithDirichletData) = f.diri_dofs

num_free_dofs(f::FESpaceWithDirichletData) = num_free_dofs(f.fespace)

num_diri_dofs(f::FESpaceWithDirichletData) = num_diri_dofs(f.fespace)

diri_tags(f::FESpaceWithDirichletData) = diri_tags(f.fespace)

function apply_constraints(
  f::FESpaceWithDirichletData, cellvec::CellVector, cellids::CellNumber)
  apply_constraints(f.fespace,cellvec,cellids)
end

function apply_constraints_rows(
  f::FESpaceWithDirichletData, cellmat::CellMatrix, cellids::CellNumber)
  apply_constraints_rows(f.fespace,cellmat,cellids)
end

function apply_constraints_cols(
  f::FESpaceWithDirichletData, cellmat::CellMatrix, cellids::CellNumber)
  apply_constraints_cols(f.fespace,cellmat,cellids)
end

function celldofids(f::FESpaceWithDirichletData)
  celldofids(f.fespace)
end

function interpolate_values(f::FESpaceWithDirichletData,fun::Function)
  free_vals, _ = interpolate_values(f.fespace,fun)
  free_vals, f.diri_dofs
end

function interpolate_values(f::FESpaceWithDirichletData{D,Z,T},val::T) where {D,Z,T}
  free_vals, _ = interpolate_values(f.fespace,val)
  free_vals, f.diri_dofs
end

function CellField(
  f::FESpaceWithDirichletData{D,Z,T},free_dofs::AbstractVector{E},diri_dofs::AbstractVector{E})where {D,Z,T,E}
  CellField(f.fespace,free_dofs,f.diri_dofs)
end

function CellBasis(f::FESpaceWithDirichletData)
  CellBasis(f.fespace)
end

Triangulation(f::FESpaceWithDirichletData) = Triangulation(f.fespace)

function interpolate_diri_values(f::FESpaceWithDirichletData, funs::Vector{<:Function})
  f.diri_dofs
end

function interpolate_diri_values(f::FESpaceWithDirichletData{D,Z,T}, vals::Vector{T}) where {D,Z,T}
  f.diri_dofs
end

end # module FESpaces
