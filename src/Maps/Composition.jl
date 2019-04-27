# Composition

function compose(f::Function,g::Field{D,S}) where {D,S}
  FieldFromCompose(f,g)
end

function compose(f::Function,g::Geomap{D,Z},u::Field{Z,S}) where {D,Z,S}
  FieldFromComposeExtended(f,g,u)
end

(∘)(f::Function,g::Field) = compose(f,g)

struct FieldFromCompose{D,O,C<:Field{D},T} <: Field{D,T}
  a::C
  op::O
end

function FieldFromCompose(f::Function,g::Field{D,S}) where {D,S}
  O = typeof(f)
  C = typeof(g)
  T = Base._return_type(f,Tuple{S})
  FieldFromCompose{D,O,C,T}(g,f)
end

function evaluate(self::FieldFromCompose{D},points::AbstractVector{Point{D}}) where D
  avals = evaluate(self.a,points)
  return broadcast(self.op,avals)
end

function gradient(self::FieldFromCompose)
  gradop = gradient(self.op)
  FieldFromCompose(gradop,self.a)
  # @santiagobadia : THIS IS WRONG
end

struct FieldFromComposeExtended{D,O,G<:Geomap{D},U<:Field,T} <: Field{D,T}
  f::O
  g::G
  u::U
end

function FieldFromComposeExtended(f::Function,g::Geomap{D,Z},u::Field{Z,S}) where {D,Z,S}
  O = typeof(f)
  G = typeof(g)
  U = typeof(u)
  T = Base._return_type(f,Tuple{Point{Z},S})
  FieldFromComposeExtended{D,O,G,U,T}(f,g,u)
end

function evaluate(self::FieldFromComposeExtended{D},points::AbstractVector{Point{D}}) where D
  gvals = evaluate(self.g,points)
  uvals = evaluate(self.u,gvals)
  return broadcast(self.f,uvals)
end

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
