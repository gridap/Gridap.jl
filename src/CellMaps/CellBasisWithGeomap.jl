attachgeomap(a::CellBasis{D},b::CellGeomap{D,D}) where D = CellBasisWithGeomap(a,b)

struct CellBasisWithGeomap{D,O,T,B<:CellBasis{D,T},G<:CellGeomap{D,D}} <: CellBasis{D,T}
  basis::B
  geomap::G
end

function CellBasisWithGeomap(basis::CellBasis{D,T},geomap::CellField{D,Point{D}},order::Int) where {D,T}
  B = typeof(basis)
  G = typeof(geomap)
  O = order
  CellBasisWithGeomap{D,O,T,B,G}(basis,geomap)
end

function CellBasisWithGeomap(basis::CellBasis{D,T},geomap::CellField{D,Point{D}}) where {D,T}
  CellBasisWithGeomap(basis,geomap,0)
end

evaluate(self::CellBasisWithGeomap{D},points::CellPoints{D}) where D = @notimplemented

gradient(self::CellBasisWithGeomap) = @notimplemented

evaluate(self::CellBasisWithGeomap{D,0},points::CellPoints{D}) where D = evaluate(self.basis,points)

function evaluate(self::CellBasisWithGeomap{D,1},points::CellPoints{D}) where D
  vals = evaluate(self.basis,points)
  jaco = gradient(self.geomap)
  jacovals = evaluate(jaco,points)
  cellnewaxis(inv(jacovals),dim=1) * vals
end

function gradient(self::CellBasisWithGeomap{D,0}) where D
  basisgrad = gradient(self.basis)
  CellBasisWithGeomap(basisgrad,self.geomap,1)
end
