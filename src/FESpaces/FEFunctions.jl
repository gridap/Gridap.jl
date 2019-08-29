module FEFunctions

using Gridap
using Gridap.Triangulations: _attach_triangulation

export FEFunction
export interpolate
import Gridap: free_dofs
import Gridap: diri_dofs
import Gridap: FESpace
import Gridap: value_type

import Gridap: evaluate, gradient, return_size
import Gridap: symmetric_gradient
import Base: div
import Gridap: trace
import Gridap: curl
import Gridap: reindex
import Gridap: restrict
import Gridap: Triangulation
import Base: size
import Base: getindex
import Base: +, -, *

import Base: zero

"""
Abstract type representing a FE Function
A FE function is a member of a FESpace.
Since a FE space is the direct sum of free + dirichlet spaces
a FE function is uniquely represented by two vectors of components aka
`free_dofs` and `diri_dofs` which are Vectors{E} with
E = eltype(T) and its fe space
"""
struct FEFunction{
  D,Z,T,E,R,
  C <: IndexCellField{Z,T}} <: IndexCellValue{R,1}
  free_dofs::AbstractVector{E}
  diri_dofs::AbstractVector{E}
  fespace::FESpace{D,Z,T}
  cellfield::C
end

free_dofs(f::FEFunction) = f.free_dofs

diri_dofs(f::FEFunction) = f.diri_dofs

FESpace(f::FEFunction) = f.fespace

Triangulation(f::FEFunction) = Triangulation(f.fespace)

value_type(::FEFunction{D,Z,T}) where {D,Z,T} = T

"""
Returns the FE function represented be the  free and dirichlet values
E = eltype(T)
"""
function FEFunction(
  fespace::FESpace,free_dofs::AbstractVector,diri_dofs::AbstractVector)
  _FEFunction(fespace,free_dofs,diri_dofs)
end

function _FEFunction(
  fespace::FESpace,free_dofs::AbstractVector,diri_dofs::AbstractVector)
  cfield = CellField(fespace,free_dofs,diri_dofs)
  FEFunction(free_dofs,diri_dofs,fespace,cfield)
end

function FEFunction(
  f::FESpaceWithDirichletData,free_vals::AbstractVector{E},diri_dofs::AbstractVector{E}) where E
  _FEFunction(f,free_vals,f.diri_dofs)
end

function FEFunction(f::FESpaceWithDirichletData,free_vals::AbstractVector)
  _FEFunction(f,free_vals,f.diri_dofs)
end

function FEFunction(
  free_dofs::AbstractVector{E},
  diri_dofs::AbstractVector{E},
  fespace::FESpace{D,Z,T},
  cellfield::CellField{Z,T}) where {D,Z,T,E}
  @assert E == eltype(T)
  R = eltype(cellfield)
  C = typeof(cellfield)
  FEFunction{D,Z,T,E,R,C}(free_dofs,diri_dofs,fespace,cellfield)
end

evaluate(f::FEFunction{D,Z},q::CellPoints{Z}) where {D,Z} = evaluate(f.cellfield,q)

for op in (:+,:-,:(gradient),:(symmetric_gradient),:(div),:(trace),:(curl))
  @eval begin
    function ($op)(a::FEFunction)
      g = $op(a.cellfield)
      trian = Triangulation(a)
      _attach_triangulation(g,trian)
    end
  end
end

for op in (:+, :-, :*)
  @eval begin

    function ($op)(a::FEFunction,b::CellField)
      trian = Triangulation(a)
      g = $op(a.cellfield,b)
      _attach_triangulation(g,trian)
    end

    function ($op)(a::CellField,b::FEFunction)
      trian = Triangulation(b)
      g = $op(a,b.cellfield)
      _attach_triangulation(g,trian)
    end

    function ($op)(a::FEFunction,b::Function)
      trian = Triangulation(a)
      cf = CellField(trian,b)
      $op(a,cf)
    end

    function ($op)(a::Function,b::FEFunction)
      trian = Triangulation(b)
      cf = CellField(trian,a)
      $op(cf,b)
    end

  end
end

return_size(f::FEFunction,s::Tuple{Int}) = return_size(f.cellfield,s)

getindex(f::FEFunction,i::Integer) = f.cellfield[i]

size(f::FEFunction) = (length(f.cellfield),)

function interpolate(this::FESpace,fun::Function)
  free_vals, diri_vals = interpolate_values(this,fun)
  FEFunction(this,free_vals,diri_vals)
end
# @santiagobadia : Using this, it always overwrites the FESpace Dirichlet
# values, if they exist. Is this what we want?

function zero(fespace::FESpace{D,Z,T}) where {D,Z,T}
  E = eltype(T)
  nf = num_free_dofs(fespace)
  nd = num_diri_dofs(fespace)
  free_vals = zeros(E,nf)
  diri_vals = zeros(E,nd)
  FEFunction(fespace,free_vals,diri_vals)
end

reindex(uh::FEFunction, indices::CellNumber{<:IndexLike}) = reindex(uh.cellfield,indices)

reindex(uh::FEFunction, indices::IndexCellNumber{<:IndexLike}) = reindex(uh.cellfield,indices)

restrict(uh::FEFunction,trian::BoundaryTriangulation) = restrict(uh.cellfield,trian)

restrict(uh::FEFunction,trian::SkeletonTriangulation) = restrict(uh.cellfield,trian)

end # module
