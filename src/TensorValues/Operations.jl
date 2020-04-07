###############################################################
# Comparison
###############################################################

(==)(a::MultiValue,b::MultiValue) = a.data == b.data
(≈)(a::MultiValue,b::MultiValue) = isapprox(collect(a.data), collect(b.data))
(≈)(a::VectorValue{0},b::VectorValue{0}) = true

function (≈)(
  a::AbstractArray{<:MultiValue}, b::AbstractArray{<:MultiValue})
  if size(a) != size(b); return false; end
  for (ai,bi) in zip(a,b)
    if !(ai≈bi); return false; end
  end
  true
end

function isless(a::VectorValue{N},b::VectorValue{N}) where N
  for d in N:-1:1
    if a[d] < b[d]
      return true
    elseif a[d] > b[d]
      return false
    else
      continue
    end
  end
  false
end

isless(a::Number,b::MultiValue) where {D,T} = all(a .< b.data)

###############################################################
# Addition / subtraction
###############################################################

for op in (:+,:-)
  @eval begin

    function ($op)(a::T) where {T<:MultiValue}
      r = map($op, a.data)
      T(r)
    end

    function ($op)(a::T where {T<:VectorValue{D}},b::T where {T<:VectorValue{D}})  where {D}
      r = broadcast(($op), a.data, b.data)
      VectorValue{D}(r)
    end

    function ($op)(a::T where {T<:TensorValue{D1,D2}},b::T where {T<:TensorValue{D1,D2}})  where {D1,D2}
      r = broadcast(($op), a.data, b.data)
      TensorValue{D1,D2}(r)
    end

  end
end

###############################################################
# Matrix Division
###############################################################

function (\)(a::T1 where {T1<:Union{TensorValue,SymTensorValue}}, b::T2) where {T2<:MultiValue}
  r = get_array(a) \ get_array(b)
  T2(r)
end

###############################################################
# Operations with other numbers
###############################################################

for op in (:+,:-,:*)
  @eval begin
    function ($op)(a::T,b::Number) where {T<:MultiValue} 
        r = broadcast($op,a.data,b)
        PT  = change_eltype(T,eltype(r))
        PT(r)
    end

    function ($op)(a::Number,b::T) where {T<:MultiValue} 
        r = broadcast($op,a,b.data)
        PT  = change_eltype(T,eltype(r))
        PT(r)
    end
  end
end

function (/)(a::T,b::Number) where {T<:MultiValue}
    r = broadcast(/,a.data,b)
    PT  = change_eltype(T,eltype(r))
    PT(r)
end

###############################################################
# Dot product (simple contraction)
###############################################################

(*)(a::VectorValue{D}, b::VectorValue{D}) where D = inner(a,b)

@generated function (*)(a::VectorValue{D1}, b::MultiValue{Tuple{D1,D2},T2,2}) where {D1,D2,T2}
    ss = String[]
    for j in 1:D2
        s = join([ "a[$i]*b[$i,$j]+" for i in 1:D1])
        push!(ss,s[1:(end-1)]*", ")
    end
    str = join(ss)
    Meta.parse("VectorValue{$D2}($str)")
end

@generated function (*)(a::MultiValue{Tuple{D1,D2},T1,2}, b::VectorValue{D1}) where {D1,D2,T1}
    ss = String[]
    for j in 1:D2
        s = join([ "b[$i]*a[$j,$i]+" for i in 1:D1])
        push!(ss,s[1:(end-1)]*", ")
    end
    str = join(ss)
    Meta.parse("VectorValue{$D2}($str)")
end

@generated function (*)(a::MultiValue{Tuple{D1,D2},T1,2}, b::MultiValue{Tuple{D3,D2},T2,2}) where {D1,D2,D3,T1,T2}
    ss = String[]
    for j in 1:D2
        for i in 1:D1
            s = join([ "a[$i,$k]*b[$k,$j]+" for k in 1:D3])
            push!(ss,s[1:(end-1)]*", ")
        end
    end
    str = join(ss)
    Meta.parse("TensorValue{$D2,$D3}(($str))")
end

@inline dot(u::VectorValue,v::VectorValue) = inner(u,v)
@inline dot(u::TensorValue,v::VectorValue) = u*v
@inline dot(u::VectorValue,v::TensorValue) = u*v

###############################################################
# Inner product (full contraction)
###############################################################

inner(a::Real,b::Real) = a*b
@generated function inner(a::MultiValue, b::MultiValue)
    @assert length(a) == length(b)
    str = join([" a[$i]*b[$i] +" for i in 1:length(a) ])
    Meta.parse(str[1:(end-1)])
end

@generated function inner(a::MultiValue{Tuple{D1,D2},T,2}, b::MultiValue{Tuple{D1,D2},T,2}) where {D1,D2,T}
    str = join([" a[$i,$j]*b[$i,$j] +" for i in 1:D1 for j in 1:D2 ])
    Meta.parse(str[1:(end-1)])
end

###############################################################
# Reductions
###############################################################

for op in (:sum,:maximum,:minimum)
    @eval begin
        $op(a::MultiValue) = $op(a.data)
    end
end

# Outer product (aka dyadic product)

"""
"""
outer(a::Real,b::Real) = a*b
outer(a::MultiValue,b::Real) = a*b
outer(a::Real,b::MultiValue) = a*b
@generated function outer(a::VectorValue{D},b::VectorValue{Z}) where {D,Z}
    str = join(["a[$i]*b[$j], " for j in 1:Z for i in 1:D])
    Meta.parse("TensorValue{$D,$Z}(($str))")
end

#@generated function outer(a::VectorValue{D},b::TensorValue{D1,D2}) where {D,D1,D2}
#  str = join(["a.array[$i]*b.array[$j,$k], "  for k in 1:D2 for j in 1:D1 for i in 1:D])
#  Meta.parse("MultiValue(SArray{Tuple{$D,$D1,$D2}}($str))")
#end

###############################################################
# Linear Algebra
###############################################################

det(a::MultiValue{Tuple{D1,D2},T,2,L}) where {D1,D2,T,L} = det(convert(SMatrix{D1,D2,T,L},a))
inv(a::MultiValue{Tuple{D1,D2},T,2,L}) where {D1,D2,T,L} = TensorValue(inv(convert(SMatrix{D1,D2,T,L},a)))

###############################################################
# Measure
###############################################################

"""
"""
meas(a::VectorValue) = sqrt(inner(a,a))
meas(a::MultiValue{S,T,2,L}) where {S,T,L} = abs(det(a))

function meas(v::TensorValue{1,2})
  n1 = v[1,2]
  n2 = -1*v[1,1]
  n = VectorValue(n1,n2)
  sqrt(n*n)
end

function meas(v::TensorValue{2,3})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1,n2,n3)
  sqrt(n*n)
end

@inline norm(u::VectorValue) = sqrt(inner(u,u))
@inline norm(u::VectorValue{0,T}) where T = sqrt(zero(T))

###############################################################
# conj
###############################################################

conj(a::VectorValue) = a
conj(a::SymTensorValue) = a
conj(a::T) where {T<:TensorValue} = T(conj(get_array(a)))

###############################################################
# Trace
###############################################################

@generated function tr(v::MultiValue{Tuple{D,D}}) where D
    str = join([" v[$i,$i] +" for i in 1:D ])
    Meta.parse(str[1:(end-1)])
end

#@generated function tr(v::MultiValue{Tuple{A,A,B}}) where {A,B}
#  str = ""
#  for k in 1:B
#    for i in 1:A
#      if i !=1
#        str *= " + "
#      end
#      str *= " v.array[$i,$i,$k]"
#    end
#    str *= ", "
#  end
#  Meta.parse("VectorValue($str)")
#end

###############################################################
# Adjoint and transpose
###############################################################

adjoint(arg::T) where {T<:TensorValue} = T(adjoint(get_array(arg)))
transpose(arg::T) where {T<:TensorValue} = T(transpose(get_array(arg)))
adjoint(arg::SymTensorValue) = arg
transpose(arg::SymTensorValue) = arg

###############################################################
# Symmetric part
###############################################################

"""
"""
@generated function symmetic_part(v::MultiValue{Tuple{D,D},T,2,L}) where {D,T,L}
    str = "("
    for j in 1:D
        for i in 1:D
            str *= "0.5*v[$i,$j] + 0.5*v[$j,$i], "
        end
    end
    str *= ")"
    Meta.parse("TensorValue($str)")
end

###############################################################
# Define new operations for Gridap types
###############################################################

for op in (:symmetic_part,)
    @eval begin
        ($op)(a::GridapType) = operate($op,a)
    end
end

for op in (:inner,:outer)
    @eval begin
        ($op)(a::GridapType,b::GridapType) = operate($op,a,b)
        ($op)(a::GridapType,b::Number)     = operate($op,a,b)
        ($op)(a::Number,    b::GridapType) = operate($op,a,b)
        ($op)(a::GridapType,b::Function)   = operate($op,a,b)
        ($op)(a::Function,  b::GridapType) = operate($op,a,b)
    end
end

