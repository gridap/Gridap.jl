
Base.:+(a::CellFunction{S,M,T,N},b::CellFunction{S,M,T,N}) where {S,M,T,N} = CellFunctionFromBaseOp(a,b,+)

Base.:-(a::CellFunction{S,M,T,N},b::CellFunction{S,M,T,N}) where {S,M,T,N} = CellFunctionFromBaseOp(a,b,-)

Base.:*(a::CellFunction{S,M,T,N},b::CellFunction{S,M,T,N}) where {S,M,T,N} = CellFunctionFromBaseOp(a,b,*)

Base.:/(a::CellFunction{S,M,T,N},b::CellFunction{S,M,T,N}) where {S,M,T,N} = CellFunctionFromBaseOp(a,b,/)

Base.:∘(f::Function,g::CellField) = compose(f,g)

function inner(a::CellField{D,T},b::CellField{D,T}) where {D,T}
  S = inner(T,T)
  CellFunctionFromInner{D,T,1,1,S,1}(a,b)
end

function inner(a::CellBasis{D,T},b::CellField{D,T}) where {D,T}
  S = inner(T,T)
  CellFunctionFromInner{D,T,2,1,S,2}(a,b)
end

function inner(a::CellBasis{D,T},b::CellBasis{D,T}) where {D,T}
  S = inner(T,T)
  CellFunctionFromInner{D,T,2,2,S,3}(a,b)
end

function expand(a::CellBasis,b::CellVector)
  CellFieldFromExpand(a,b)
end

# @fverdugo TODO: optimize for ConstantCellArrays
function compose(f::Function,g::CellField{D,S}) where {D,S}
  O = typeof(f)
  C = typeof(g)
  T = f(S)
  CellFieldFromCompose{D,O,C,T}(g,f)
end

inner(a::CellFieldValues{T},b::CellFieldValues{T}) where T = binner(a,b)

inner(a::CellBasisValues{T},b::CellFieldValues{T}) where T = binner(a,cellnewaxis(b,dim=1))

inner(a::CellBasisValues{T},b::CellBasisValues{T}) where T = binner(cellnewaxis(a,dim=2),cellnewaxis(b,dim=1))

expand(a::CellBasisValues,b::CellFieldValues) = cellsum(bouter(a,cellnewaxis(b,dim=2)),dim=1)

# Ancillary types associated with operations above

struct CellFieldFromCompose{D,O,C<:CellField{D},T} <: CellField{D,T}
  a::C
  op::O
end

function evaluate(self::CellFieldFromCompose{D,O,C,T},points::CellPoints{D}) where {D,O,C,T}
  avals = evaluate(self.a,points)
  CellArrayFromGivenUnaryOp{O,typeof(avals),T,1}(avals,self.op)
end

function gradient(self::CellFieldFromCompose{D,O,C,T}) where {D,O,C,T}
  gradop = gradient(self.op)
  OG = typeof(gradop)
  TG = gradop(Point{D})
  CellFieldFromCompose{D,OG,C,TG}(self.a,gradop)
end

# @fverdugo this is type unstable (however, not a big problem here)
struct CellFunctionFromBaseOp{O,S,M,T,N} <: CellFunction{S,M,T,N}
  a::CellFunction{S,M,T,N}
  b::CellFunction{S,M,T,N}
  op::O
end

function evaluate(self::CellFunctionFromBaseOp{O,S,M,T,N},input::CellArray{S,M}) where {O,S,M,T,N}
  avals = evaluate(self.a,input)
  bvals = evaluate(self.b,input)
  self.op(avals,bvals)
end

struct CellFunctionFromInner{D,R,A,B,T,N} <: CellFunction{Point{D},1,T,N}
  a::CellFunction{Point{D},1,R,A}
  b::CellFunction{Point{D},1,R,B}
end

function evaluate(self::CellFunctionFromInner{D},points::CellPoints{D}) where D
  a = evaluate(self.a,points)
  b = evaluate(self.b,points)
  inner(a,b)
end

struct CellFieldFromExpand{D,S,R,T} <: CellField{D,T}
  basis::CellBasis{D,S}
  coeffs::CellVector{R}
end

function CellFieldFromExpand(basis::CellBasis{D,S},coeffs::CellVector{R}) where {D,S,R}
  T = outer(S,R)
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

