
# Unary operations

for op in (:+, :-)
  @eval begin
    function ($op)(a::CellFunction)
      CellFunctionFromUnaryOp($op,a)
    end
  end
end

struct CellFunctionFromUnaryOp{O,A,S,M,T,N} <: CellFunction{S,M,T,N}
  op::O
  a::A
end

function CellFunctionFromUnaryOp(op::Function,a::CellFunction{S,M,T,N}) where {S,M,T,N}
  O = typeof(op)
  A = typeof(a)
  R = Base._return_type(op,Tuple{T})
  CellFunctionFromUnaryOp{O,A,S,M,R,N}(op,a)
end

function evaluate(self::CellFunctionFromUnaryOp{O,A,S,M,T,N},input::CellArray{S,M}) where {O,A,S,M,T,N}
  avals = evaluate(self.a,input)
  self.op(avals)
end

# Binary operations

for op in (:+, :-, :*, :/)
  @eval begin
    function ($op)(a::CellFunction,b::CellFunction)
      CellFunctionFromBinaryOp($op,a,b)
    end
  end
end

function inner(a::CellField{D,T},b::CellField{D,T}) where {D,T}
  CellFunctionFromBinaryOp(varinner,a,b)
end

function inner(a::CellBasis{D,T},b::CellField{D,T}) where {D,T}
  CellFunctionFromBinaryOp(varinner,a,b)
end

function inner(a::CellBasis{D,T},b::CellBasis{D,T}) where {D,T}
  CellFunctionFromBinaryOp(varinner,a,b)
end

function expand(a::CellBasis,b::CellVector)
  CellFieldFromExpand(a,b)
end

struct CellFunctionFromBinaryOp{O,A,B,S,M,T,N} <: CellFunction{S,M,T,N}
  op::O
  a::A
  b::B
end

function CellFunctionFromBinaryOp(op::Function,a::CellFunction{S,M,TA,NA},b::CellFunction{S,M,TB,NB}) where {S,M,TA,NA,TB,NB}
  O = typeof(op)
  A = typeof(a)
  B = typeof(b)
  R = Base._return_type(op,Tuple{TA,TB})
  N = max(NA,NB)
  CellFunctionFromBinaryOp{O,A,B,S,M,R,N}(op,a,b)
end

function evaluate(self::CellFunctionFromBinaryOp{O,A,B,S,M,T,N},input::CellArray{S,M}) where {O,A,B,S,M,T,N}
  avals = evaluate(self.a,input)
  bvals = evaluate(self.b,input)
  self.op(avals,bvals)
end

struct CellFieldFromExpand{D,S,R,T} <: CellField{D,T}
  basis::CellBasis{D,S}
  coeffs::CellVector{R}
end

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

struct CellFieldFromCompose{D,O,C<:CellField{D},T} <: CellField{D,T}
  a::C
  op::O
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
  ConstantCellArray(b,length(avals))
end

function gradient(self::CellFieldFromCompose)
  gradop = gradient(self.op)
  CellFieldFromCompose(gradop,self.a)
end

struct CellFieldFromComposeExtended{D,O,G<:CellGeomap{D},U<:CellField,T} <: CellField{D,T}
  f::O
  g::G
  u::U
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

