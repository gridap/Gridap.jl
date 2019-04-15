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
