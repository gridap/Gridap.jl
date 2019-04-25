# Unary operations

for op in (:+, :-)
  @eval begin
    function ($op)(a::CellMap)
      CellMapFromUnaryOp($op,a)
    end
  end
end

struct CellMapFromUnaryOp{O,A,S,M,T,N} <: IterCellMap{S,M,T,N}
  op::O
  a::A
end

@inline function Base.iterate(this::CellMapFromUnaryOp{O,A,S,M,T,N}) where {O,A,S,M,T,N}
  anext = iterate(this.a)
  iteratekernel(this,anext)
end

@inline function Base.iterate(this::CellMapFromUnaryOp{O,A,S,M,T,N},state) where {O,A,S,M,T,N}
  v, astate = state
  anext = iterate(this.a,astate)
  iteratekernel(this,anext)
end

function iteratekernel(this::CellMapFromUnaryOp,anext)
  if anext === nothing; return nothing end
  a, astate = anext
  v = MapFromUnaryOp(this.op,a)
  state = (v, astate)
  (v, state)
end

function CellMapFromUnaryOp(op::Function,a::CellMap{S,M,T,N}) where {S,M,T,N}
  O = typeof(op)
  A = typeof(a)
  R = Base._return_type(op,Tuple{T})
  CellMapFromUnaryOp{O,A,S,M,R,N}(op,a)
end

function evaluate(self::CellMapFromUnaryOp{O,A,S,M,T,N},input::CellArray{S,M}) where {O,A,S,M,T,N}
  avals = evaluate(self.a,input)
  self.op(avals)
end

# Binary operations

for op in (:+, :-, :*, :/)
  @eval begin
    function ($op)(a::CellMap,b::CellMap)
      CellMapFromBinaryOp($op,a,b)
    end
  end
end

function inner(a::CellField{D,T},b::CellField{D,T}) where {D,T}
  CellMapFromBinaryOp(varinner,a,b)
end

function inner(a::CellBasis{D,T},b::CellField{D,T}) where {D,T}
  CellMapFromBinaryOp(varinner,a,b)
end

function inner(a::CellBasis{D,T},b::CellBasis{D,T}) where {D,T}
  CellMapFromBinaryOp(varinner,a,b)
end

function expand(a::CellBasis,b::CellVector)
  CellFieldFromExpand(a,b)
end

struct CellMapFromBinaryOp{O,A,B,S,M,T,N} <: IterCellMap{S,M,T,N}
  op::O
  a::A
  b::B
end

@inline function Base.iterate(this::CellMapFromBinaryOp{O,A,B,S,M,T,N}) where {O,A,B,S,M,T,N}
  anext = iterate(this.a)
  bnext = iterate(this.b)
  iteratekernel(this,anext,bnext)
end

@inline function Base.iterate(this::CellMapFromBinaryOp{O,A,B,S,M,T,N},state) where {O,A,B,S,M,T,N}
  v, astate, bstate = state
  anext = iterate(this.a,astate)
  bnext = iterate(this.b,bstate)
  iteratekernel(this,anext,bnext)
end

function iteratekernel(this::CellMapFromBinaryOp,anext,bnext)
  if anext === nothing; return nothing end
  if bnext === nothing; return nothing end
  a, astate = anext
  b, bstate = bnext
  v = MapFromBinaryOp(this.op,a,b)
  state = (v, astate, bstate)
  (v, state)
end

function CellMapFromBinaryOp(op::Function,a::CellMap{S,M,TA,NA},b::CellMap{S,M,TB,NB}) where {S,M,TA,NA,TB,NB}
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{TA,TB})
  N = max(NA,NB)
  CellMapFromBinaryOp{O,A,B,S,M,R,N}(op,a,b)
end

function evaluate(self::CellMapFromBinaryOp{O,A,B,S,M,T,N},input::CellArray{S,M}) where {O,A,B,S,M,T,N}
  avals = evaluate(self.a,input)
  bvals = evaluate(self.b,input)
  self.op(avals,bvals)
end

struct CellFieldFromExpand{D,S,R,T<:FieldValue} <: IterCellField{D,T}
  basis::CellBasis{D,S}
  coeffs::CellVector{R}
end


@inline function Base.iterate(this::CellFieldFromExpand{D,S,R,T}) where {D,S,R,T}
  bnext = iterate(this.basis)
  cnext = iterate(this.coeffs)
  iteratekernel(this,bnext,cnext)
end

@inline function Base.iterate(this::CellFieldFromExpand{D,S,R,T},state) where {D,S,R,T}
  v, bstate, cstate = state
  bnext = iterate(this.basis,bstate)
  cnext = iterate(this.coeffs,cstate)
  iteratekernel(this,bnext,cnext)
end

function iteratekernel(this::CellFieldFromExpand,bnext,cnext)
  if bnext === nothing; return nothing end
  if cnext === nothing; return nothing end
  b, bstate = bnext
  c, cstate = cnext
  v = FieldFromExpand(b,c)
  state = (v, bstate, cstate)
  (v, state)
end

# @santiagobadia : NOT IMPLEMENTED IN THE OTHER STRUCTS!!! I want to talk to
# @fverdugo first

function CellFieldFromExpand(basis::CellBasis{D,S},coeffs::CellVector{R}) where {D,S,R}
  T = Base._return_type(outer,Tuple{S,R})
  CellFieldFromExpand{D,S,R,T}(basis,coeffs)
end

function evaluate(self::CellFieldFromExpand{D},points::CellPoints{D}) where D
  basisvals = evaluate(self.basis,points)
  expand(basisvals,self.coeffs)
end

function gradient(self::CellFieldFromExpand)
  gradbasis = gradient(self.basis)
  CellFieldFromExpand(gradbasis,self.coeffs)
end

varinner(a::FieldValue,b::FieldValue) = inner(a,b)

varinner(a::CellFieldValues{T},b::CellFieldValues{T}) where T = inner(a,b)

varinner(a::CellBasisValues{T},b::CellFieldValues{T}) where T = inner(a,cellnewaxis(b,dim=1))

varinner(a::CellBasisValues{T},b::CellBasisValues{T}) where T = inner(cellnewaxis(a,dim=2),cellnewaxis(b,dim=1))

expand(a::CellBasisValues,b::CellFieldValues) = cellsum(outer(a,cellnewaxis(b,dim=2)),dim=1)

# Composition

function compose(f::Function,g::CellField{D,S}) where {D,S}
  CellFieldFromCompose(f,g)
end

function compose(f::Function,g::CellGeomap{D,Z},u::CellField{Z,S}) where {D,Z,S}
  CellFieldFromComposeExtended(f,g,u)
end

(∘)(f::Function,g::CellField) = compose(f,g)

struct CellFieldFromCompose{D,O,C<:CellField{D},T} <: IterCellField{D,T}
  a::C
  op::O
end

@inline function Base.iterate(this::CellFieldFromCompose{D,O,C,T}) where {D,O,C,T}
  anext = iterate(this.a)
  iteratekernel(this,anext)
end

@inline function Base.iterate(this::CellFieldFromCompose{D,O,C,T},state) where {D,O,C,T}
  v, astate = state
  anext = iterate(this.a,astate)
  iteratekernel(this,anext)
end

function iteratekernel(this::CellFieldFromCompose,anext)
  if anext === nothing; return nothing end
  a, astate = anext
  v = FieldFromCompose(this.op,a)
  state = (v, astate)
  (v, state)
end

function CellFieldFromCompose(f::Function,g::CellField{D,S}) where {D,S}
  O = typeof(f)
  C = typeof(g)
  T = Base._return_type(f,Tuple{S})
  CellFieldFromCompose{D,O,C,T}(g,f)
end

function evaluate(self::CellFieldFromCompose{D},points::CellPoints{D}) where D
  avals = evaluate(self.a,points)
  composekernel(self,avals)
end

function composekernel(self::CellFieldFromCompose,avals::CellArray)
  CellArrayFromBroadcastUnaryOp(self.op,avals)
end

function composekernel(self::CellFieldFromCompose,avals::ConstantCellArray)
  b = broadcast(self.op,celldata(avals))
  ConstantCellValue(b,length(avals))
end

function gradient(self::CellFieldFromCompose)
  gradop = gradient(self.op)
  CellFieldFromCompose(gradop,self.a)
end

struct CellFieldFromComposeExtended{D,O,G<:CellGeomap{D},U<:CellField,T} <: IterCellField{D,T}
  f::O
  g::G
  u::U
end

@inline function Base.iterate(this::CellFieldFromComposeExtended{D,O,G,U,T}) where {D,O,G,U,T}
  gnext = iterate(this.g)
  unext = iterate(this.u)
  iteratekernel(this,gnext,unext)
end

@inline function Base.iterate(this::CellFieldFromComposeExtended{D,O,G,U,T},state) where {D,O,G,U,T}
  v, gstate, ustate = state
  gnext = iterate(this.g,astate)
  unext = iterate(this.u,bstate)
  iteratekernel(this,gnext,unext)
end

function iteratekernel(this::CellFieldFromComposeExtended,gnext,unext)
  if gnext === nothing; return nothing end
  if unext === nothing; return nothing end
  g, gstate = gnext
  u, ustate = unext
  v = FieldFromComposeExtended(this.op,g,u)
  state = (v, gstate, ustate)
  (v, state)
end

function CellFieldFromComposeExtended(f::Function,g::CellGeomap{D,Z},u::CellField{Z,S}) where {D,Z,S}
  O = typeof(f)
  G = typeof(g)
  U = typeof(u)
  T = Base._return_type(f,Tuple{Point{Z},S})
  CellFieldFromComposeExtended{D,O,G,U,T}(f,g,u)
end

function evaluate(self::CellFieldFromComposeExtended{D},points::CellPoints{D}) where D
  gvals = evaluate(self.g,points)
  uvals = evaluate(self.u,gvals)
  CellArrayFromBroadcastBinaryOp(self.f,gvals,uvals)
end

function gradient(self::CellFieldFromComposeExtended)
  gradf = gradient(self.f)
  CellFieldFromComposeExtended(gradf,self.g,self.u)
end
