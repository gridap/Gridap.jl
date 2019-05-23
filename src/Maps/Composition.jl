# Composition

(∘)(f::Function,g::Field) = compose(f,g)

function compose(f::Function,g::Field{D,S}) where {D,S}
  FieldFromCompose(f,g)
end

function compose(f::Function,g::Geomap{D,Z},u::Field{Z,S}) where {D,Z,S}
  FieldFromComposeExtended(f,g,u)
end

function attachgeomap(a::Basis{D,T},b::Geomap{D,D}) where {D,T}
  BasisWithGeoMap(a,b)
end

# Helpers

struct BasisWithGeoMap{
  O,D,T,
  DD,A<:Basis{D,T},B<:Field{D,TensorValue{D,DD}}} <: MapFromBinaryOp{Point{D},1,T,2}
  map::A
  jaco::B
  acache::CachedMatrix{T,Matrix{T}}
  bcache::CachedVector{TensorValue{D,DD},Vector{TensorValue{D,DD}}}
end

function BasisWithGeoMap(
  a::Basis{D,T},
  b::Geomap{D,D}) where {D,T}
  jaco = gradient(b)
  BasisWithGeoMap(a,jaco,0)
end

function BasisWithGeoMap(
  a::Basis{D,T},
  jaco::Field{D,<:TensorValue{D}},
  order::Int) where {D,T}

  O = order
  DD = D*D
  A = typeof(a)
  B = typeof(jaco)
  acache = CachedMatrix(T)
  bcache = CachedVector(TensorValue{D,DD})
  BasisWithGeoMap{O,D,T,DD,A,B}(a,jaco,acache,bcache)
end

leftmap(self::BasisWithGeoMap) = self.map

rightmap(self::BasisWithGeoMap) = self.jaco

computesize(self::BasisWithGeoMap, sa, sb) = sa

function computevals!(self::BasisWithGeoMap, a, b, v)
  ndofs, npoints = size(a)
  for j in 1:npoints
    for i in 1:ndofs
    v[i,j] = inv(b[j])*a[i,j]
    end
  end
end

leftcachedarray(self::BasisWithGeoMap) = self.acache

rightcachedarray(self::BasisWithGeoMap) = self.bcache

gradient(self::BasisWithGeoMap) = @notimplemented

function gradient(self::BasisWithGeoMap{0})
  g = gradient(self.map)
  BasisWithGeoMap(g,self.jaco,1)
end

function evaluate!(
  self::BasisWithGeoMap{0,D,T},
  points::AbstractVector{Point{D}},
  v::AbstractMatrix{T}) where {D,T}
  evaluate!(self.map,points,v)
end

const FieldFromCompose{D,O,T,U,A} =  MapFromBroadcastUnaryOp{O,Point{D},1,T,1,U,A}

function FieldFromCompose(f::Function,g::Field{D,S}) where {D,S}
  O = typeof(f)
  A = typeof(g)
  T = Base._return_type(f,Tuple{S})
  cache = CachedArray(S,1)
  FieldFromCompose{D,O,T,S,A}(f,g,cache)
end

function gradient(self::FieldFromCompose)
  gradop = gradient(self.op)
  FieldFromCompose(gradop,inputmap(self))
end

struct FieldFromComposeExtended{
  D,T,O,Z,S,G<:Geomap{D,Z},U<:Field{Z,S}} <: Field{D,T}
  f::O
  g::G
  u::U
  gcache::CachedVector{Point{Z},Vector{Point{Z}}}
  ucache::CachedVector{S,Vector{S}}
end

function FieldFromComposeExtended(
  f::Function,g::Geomap{D,Z},u::Field{Z,S}) where {D,Z,S}
  O = typeof(f)
  G = typeof(g)
  U = typeof(u)
  T = Base._return_type(f,Tuple{Point{Z},S})
  gcache = CachedArray(Point{Z},1)
  ucache = CachedArray(S,1)
  FieldFromComposeExtended{D,T,O,Z,S,G,U}(f,g,u,gcache,ucache)
end

function evaluate!(
  self::FieldFromComposeExtended{D,T},
  points::AbstractVector{Point{D}},
  v::AbstractVector{T}) where {D,T}
  setsize!(self.gcache,size(points))
  setsize!(self.ucache,size(points))
  evaluate!(self.g,points,self.gcache)
  evaluate!(self.u,points,self.ucache)
  broadcast!(self.f,v,self.gcache,self.ucache)
end

return_size(self::FieldFromComposeExtended,s::Tuple{Int}) = s

function gradient(self::FieldFromComposeExtended)
  gradf = gradient(self.f)
  FieldFromComposeExtended(gradf,self.g,self.u)
end

# @santiagobadia : I think the following struct would be useful

# struct MapComposition{A<:Map{Q,O,T,N} where {Q,O}, B<:Map{S,M,Q,O} where {Q,O}} <: Map{S,M,T,N}
#   map_left::A # A∘B(x)
#   map_right::B
# end
#
# function compose(f::Map,g::Map)
#   MapComposition(f,g)
# end
#
# (∘)(f::Map,g::Map) = compose(f,g)
#
# function evaluate(this::MapComposition{S,M,T,N},points::AbstractArray{S,M}) where {S,M,T,N}
#   bvals = evaluate(self.map_right,points)
#   avals = evaluate(self.map_left,gvals)
# end
#
# function gradient(this::MapComposition)
#   @notimplemented
#   # Here, the chain rule
# end
