module CellBasesWithGeomap

using Gridap
using Gridap.Helpers
using Gridap.Maps
using Gridap.FieldValues
using Gridap.Maps: BasisWithGeoMap
using Gridap.CellMaps
using Gridap.CellValues

import Gridap.Maps: attachgeomap
import Base: iterate, length
import Gridap: evaluate, gradient
import Gridap: return_size

function attachgeomap(a::CellBasis{D},b::CellGeomap{D,D}) where D
  CellBasisWithGeomap(a,b)
end

struct CellBasisWithGeomap{
  O,D,T,B<:CellBasis{D,T},J<:CellField{D,<:TensorValue{D}},R<:BasisWithGeoMap} <: IterCellBasis{D,T,R}
  basis::B
  jaco::J
end

function CellBasisWithGeomap(
  basis::CellBasis{D}, geomap::CellGeomap{D,D}) where D
  jaco = gradient(geomap)
  CellBasisWithGeomap(basis,jaco,0)
end

function CellBasisWithGeomap(
  basis::CellBasis{D,T},
  jaco::CellField{D,<:TensorValue{D}},
  order::Int) where {D,T}
  @assert length(basis) == length(jaco)
  B = typeof(basis)
  J = typeof(jaco)
  O = order
  DD = D*D
  R = BasisWithGeoMap{O,D,T,DD,eltype(basis),eltype(jaco)}
  CellBasisWithGeomap{O,D,T,B,J,R}(basis,jaco)
end

length(self::CellBasisWithGeomap) = length(self.basis)

@inline function Base.iterate(this::CellBasisWithGeomap)
  jnext = iterate(this.jaco)
  bnext = iterate(this.basis)
  iteratekernel(this,jnext,bnext)
end

@inline function Base.iterate(this::CellBasisWithGeomap, state)
  v, jstate, bstate = state
  jnext = iterate(this.jaco,jstate)
  bnext = iterate(this.basis,bstate)
  iteratekernel(this,jnext,bnext)
end

function iteratekernel(this::CellBasisWithGeomap{O},jnext,bnext) where O
  if jnext === nothing; return nothing end
  if bnext === nothing; return nothing end
  j, jstate = jnext
  b, bstate = bnext
  v = BasisWithGeoMap(b,j,O) # TODO this allocates a new object
  state = (v, jstate, bstate)
  (v, state)
end

function evaluate(self::CellBasisWithGeomap{O,D},points::CellPoints{D}) where {O,D}
  @notimplemented
end

function evaluate(self::CellBasisWithGeomap{0,D},points::CellPoints{D}) where D
  evaluate(self.basis,points)
end

function evaluate(self::CellBasisWithGeomap{1,D},points::CellPoints{D}) where D
  vals = evaluate(self.basis,points)
  jacovals = evaluate(self.jaco,points)
  cellnewaxis(inv(jacovals),dim=1) * vals
end

gradient(self::CellBasisWithGeomap) = @notimplemented

function gradient(self::CellBasisWithGeomap{0})
  basisgrad = gradient(self.basis)
  CellBasisWithGeomap(basisgrad,self.jaco,1)
end

return_size(self::CellBasisWithGeomap,s::Tuple{Int}) = return_size(self.basis,s)

end # module CellBasesWithGeomap
