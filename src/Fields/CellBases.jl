"""
Concrete implementation of CellBasis for the case of the same interpolation on all cells
"""
struct CellBasisFromSingleInterpolation{D,T,B<:Basis{D,T}} <: CellBasis{D,T}
  basis::B
  function CellBasisFromSingleInterpolation{D,T,B}(basis::Basis{D,T}) where {D,T,B}
    new(basis)
  end
end

function CellBasisFromSingleInterpolation(basis::Basis{D,T}) where {D,T}
  CellBasisFromSingleInterpolation{D,T,typeof(basis)}(basis)
end

function evaluate(self::CellBasisFromSingleInterpolation{D,T},points::ConstantCellArray{Point{D},1}) where {D,T}
  ndofs = length(self.basis)
  npoints = length(points.array)
  values = Array{T,2}(undef,(ndofs,npoints))
  evaluate!(self.basis,points.array,values)
  ConstantCellArray(values,points.length)
end

function evaluate(self::CellBasisFromSingleInterpolation{D},points::CellPoints{D}) where D
  CellBasisValuesFromSingleInterpolation(self.basis,points)
end

function gradient(self::CellBasisFromSingleInterpolation)
  gradbasis = gradient(self.basis)
  CellBasisFromSingleInterpolation(gradbasis)
end

# Ancillary types

# @fverdugo a way to say the following?
# B <: Basis{D,T}
# C <: CellPoints{D}
# T <: FieldValue
struct CellBasisValuesFromSingleInterpolation{D,T,B,C} <: CellArrayFromUnaryOp{C,T,2}
  basis::B
  ndofs::Int
  cellpoints::C
end

function CellBasisValuesFromSingleInterpolation(
  basis::Basis{D,T},cellpoints::CellPoints{D}) where {D,T}
  B = typeof(basis)
  C = typeof(cellpoints)
  CellBasisValuesFromSingleInterpolation{D,T,B,C}(basis,length(basis),cellpoints)
end

inputcellarray(self::CellBasisValuesFromSingleInterpolation) = self.cellpoints

function computesize(self::CellBasisValuesFromSingleInterpolation, asize)
  (self.ndofs,asize[1])
end

function computevals!(self::CellBasisValuesFromSingleInterpolation, a, v)
  evaluate!(self.basis,a,v)
end
