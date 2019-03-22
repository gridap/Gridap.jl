# Concrete implementations

"""
Concrete implementation for the case of the same interpolation on
all cells, but arbitrary sampling points in each cell. This is typically needed
for unfitted methods
"""
struct CellBasisValuesFromSingleInterpolation{T,D} <: CellBasisValues{T}
  basis::MultivariatePolynomialBasis{D,T}
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

maxsize(self::CellBasisValuesFromSingleInterpolation) = (length(self.basis), maxlength(self.points))

"""
Concrete implementation for the case of the same interpolation
and the same sampling points on all cells
"""
struct ConstantCellBasisValues{D,T} <: IndexableCellArray{T,2}
  basis::MultivariatePolynomialBasis{D,T}
  points::Array{Point{D},1}
  l::Int
  values::Array{T,2}
end

function ConstantCellBasisValues(
  basis::MultivariatePolynomialBasis{D,T}, points::Array{Point{D},1}, l::Int) where {D,T}
  ndofs = length(basis)
  npoin = length(points)
  values = Array{T,2}(undef,(ndofs,npoin))
  evaluate!(basis,points,values)
  ConstantCellBasisValues(basis,points,l,values)
end

Base.length(self::ConstantCellBasisValues) = self.l

maxsize(self::ConstantCellBasisValues) = (length(self.basis),length(self.points))

Base.getindex(self::ConstantCellBasisValues,cell::Int) = self.values


struct CellBasisFromSingleInterpolation{D,T} <: CellBasis{D,T}
  basis::MultivariatePolynomialBasis{D,T}
end

function evaluate(
  self::CellBasisFromSingleInterpolation{D,T}, cellpoints::CellPoints{D}) where {D,T}
  if isa(cellpoints,ConstantCellArray)
    points = cellpoints.array
    l = length(cellpoints)
    ConstantCellBasisValues(self.basis,points,l)
  else
    CellBasisValuesFromSingleInterpolation(self.basis,cellpoints)
  end
end

function gradient(self::CellBasisFromSingleInterpolation)
  grad_basis = gradient(self.basis)
  CellBasisFromSingleInterpolation(grad_basis)
end

"""
Values of the gradients of a cell basis mapped with a given Jacobian
J stands for the type of the values of the Jacobian, i.e., gradient(Point{D},Val(D))
"""
struct GradOfCellBasisValuesMappedWithJacob{T,J} <: CellBasisValues{T}
  gradrefbasis::CellBasisValues{T}
  jacobian::CellFieldValues{J}
end

function Base.length(self::GradOfCellBasisValuesMappedWithJacob)
  @assert length(self.gradrefbasis) == length(self.jacobian)
  length(self.jacobian)
end

maxsize(self::GradOfCellBasisValuesMappedWithJacob) = maxsize(self.gradrefbasis)

function Base.iterate(self::GradOfCellBasisValuesMappedWithJacob{T,J}) where {T,J}
  (maxdofs,maxqpoints) = maxsize(self.gradrefbasis)
  values = Array{T,2}(undef,(maxdofs,maxqpoints))
  basisvalsnext = iterate(self.gradrefbasis)
  jacobianvalsnext = iterate(self.jacobian)
  state = (values,basisvalsnext,jacobianvalsnext)
  iterate(self,state)
end

function Base.iterate(self::GradOfCellBasisValuesMappedWithJacob,state)
  (values,basisvalsnext,jacobianvalsnext) = state
  if basisvalsnext == nothing || jacobianvalsnext == nothing
    nothing
  else
    (basisvals,basisvalsstate) = basisvalsnext
    (jacobianvals,jacobianvalsstate) = jacobianvalsnext
    mult_by_inverse_kernel!(basisvals,jacobianvals,values)
    npoints = size(basisvals,2)
    ndofs = size(basisvals,1)
    vals = @view values[1:ndofs,1:npoints]
    basisvalsnext = iterate(self.gradrefbasis,basisvalsstate)
    jacobianvalsnext = iterate(self.jacobian,jacobianvalsstate)
    state = (values,basisvalsnext,jacobianvalsnext)
    (vals, state)
  end
end

function mult_by_inverse_kernel!(basisvals,jacobianvals,mappedbasisvals)
  @assert size(basisvals,2) == length(jacobianvals)
  @assert size(basisvals,1) <= size(mappedbasisvals,1)
  @assert size(basisvals,2) <= size(mappedbasisvals,2)
  ndofs = size(basisvals,1)
  npoin = size(basisvals,2)
  @inbounds for i = 1:npoin
    invJ = inv(jacobianvals[i])
    @inbounds for j = 1:ndofs
      mappedbasisvals[j,i] = invJ*basisvals[j,i]
    end
  end
end

"""
Un-evaluated version of GradOfCellBasisValuesMappedWithJacob
"""
struct GradOfCellBasisMappedWithJacob{D,T,J} <: CellBasis{D,T}
  gradrefbasis::CellBasis{D,T}
  jacobian::CellField{D,J}
end

function evaluate(
  self::GradOfCellBasisMappedWithJacob{D,T,J},points::CellPoints{D}) where {D,T,J}
  basisvals = evaluate(self.gradrefbasis,points)
  jacobianvals = evaluate(self.jacobian,points)
  GradOfCellBasisValuesMappedWithJacob(basisvals,jacobianvals)
end

gradient(::GradOfCellBasisMappedWithJacob) = @notimplemented

"""
This type implements the result of mapderivatives
"""
struct CellBasisWithMappedDerivatives{D,T} <: CellBasis{D,T}
  basis::CellBasis{D,T}
  geomap::CellField{D,Point{D}}
end

evaluate(self::CellBasisWithMappedDerivatives,points) = evaluate(self.basis,points)

function gradient(self::CellBasisWithMappedDerivatives{D,T}) where {D,T}
  gradrefbasis = gradient(self.basis)
  jacobian = gradient(self.geomap)
  TG = gradient(T,Val(D))
  J = gradient(Point{D},Val(D))
  GradOfCellBasisMappedWithJacob{D,TG,J}(gradrefbasis,jacobian)
end
