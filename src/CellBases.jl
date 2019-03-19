export CellBasis
export evaluate

abstract type CellBasis{T} end

evaluate(::CellBasis{T} where T,::CellPoints{D} where D)::CellBasisValues{T} = @abstractmethod

# Concrete implementations

"""
Concrete implementation for the case of the same interpolation on
all cells, but arbitrary sampling points in each cell. This is typically needed
for unfitted methods
"""
struct CellBasisValuesFromSingleInterpolation{T,D} <: CellBasisValues{T}
  basis::MultivariatePolynomialBasis{T,D}
  points::CellPoints{D}
end

Base.length(self::CellBasisValuesFromSingleInterpolation) = length(self.points)

function Base.iterate(self::CellBasisValuesFromSingleInterpolation{T,D}) where {T,D}
  maxpoints = maxlength(self.points)
  ndofs = length(self.basis)
  values = Array{T,2}(undef,(ndofs,maxpoints))
  pnext = iterate(self.points)
  state = (values,pnext)
  iterate(self,state)
end

function Base.iterate(self::CellBasisValuesFromSingleInterpolation,state)
  (values,pnext) = state
  if pnext == nothing
    nothing
  else
    (points,pstate) = pnext
    evaluate!(self.basis,points,values)
    vals = @view values[:,1:length(points)]
    pnext = iterate(self.points,pstate)
    state = (values,pnext)
    (vals, state)
  end
end

maxlength(self::CellBasisValuesFromSingleInterpolation) = maxlength(self.points)*length(self.basis)

"""
Concrete implementation for the case of the same interpolation
and the same sampling points on all cells
"""
struct ConstantCellBasisValues{T,D} <: IndexableCellArray{T,2}
  basis::MultivariatePolynomialBasis{T,D}
  points::Array{Point{D},1}
  l::Int
  values::Array{T,2}
end

function ConstantCellBasisValues(
  basis::MultivariatePolynomialBasis{T,D},
  points::Array{Point{D},1},
  l::Int) where {T,D}
  ndofs = length(basis)
  npoin = length(points)
  values = Array{T,2}(undef,(ndofs,npoin))
  evaluate!(basis,points,values)
  ConstantCellBasisValues(basis,points,l,values)
end

Base.length(self::ConstantCellBasisValues) = self.l

maxlength(self::ConstantCellBasisValues) = length(self.points)*length(self.basis)

Base.getindex(self::ConstantCellBasisValues,cell::Int) = self.values

struct CellBasisFromSingleInterpolation{T,D} <: CellBasis{T}
  basis::MultivariatePolynomialBasis{T,D}
end

function evaluate(
  self::CellBasisFromSingleInterpolation{T,D},
  cellpoints::CellPoints{D}) where {T,D}
  if isa(cellpoints,ConstantCellArray)
    points = cellpoints[1]
    l = length(cellpoints)
    ConstantCellBasisValues(self.basis,points,l)
  else
    CellBasisValuesFromSingleInterpolation(self.basis,cellpoints)
  end
end

