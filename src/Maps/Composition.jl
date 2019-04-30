# Composition

function compose(f::Function,g::Field{D,S}) where {D,S}
  FieldFromCompose(f,g)
end

function compose(f::Function,g::Geomap{D,Z},u::Field{Z,S}) where {D,Z,S}
  FieldFromComposeExtended(f,g,u)
end

(∘)(f::Function,g::Field) = compose(f,g)

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

#struct FieldFromCompose{D,O,C<:Field{D},T} <: Field{D,T}
#  a::C
#  op::O
#  cache::CachedVector{T,Vector{T}}
#end
#
#function FieldFromCompose(f::Function,g::Field{D,S}) where {D,S}
#  O = typeof(f)
#  C = typeof(g)
#  T = Base._return_type(f,Tuple{S})
#  cache = CachedArray(S,1)
#  FieldFromCompose{D,O,C,T}(g,f,cache)
#end
#
#function evaluate!(
#  self::FieldFromCompose{D,O,C,T},
#  points::AbstractVector{Point{D}},
#  v::AbstractVector{T}) where {D,O,C,T}
#  setsize!(self.cache,size(points))
#  evaluate!(self.a,points,self.cache)
#  broadcast!(self.op,v,self.cache)
#end
#
#return_size(self::FieldFromCompose, s::NTuple{1,Int}) = s
#
#function gradient(self::FieldFromCompose)
#  gradop = gradient(self.op)
#  FieldFromCompose(gradop,self.a)
#  # @santiagobadia : THIS IS WRONG
#end

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
  evaluate!(self.u,self.gcache,self.ucache)
  broadcast!(self.f,v,self.gcache,self.ucache)
end

return_size(self::FieldFromComposeExtended,s::Tuple{Int}) = s

function gradient(self::FieldFromComposeExtended)
  gradf = gradient(self.f)
  FieldFromComposeExtended(gradf,self.g,self.u)
  # @santiagobadia : THIS IS WRONG
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
