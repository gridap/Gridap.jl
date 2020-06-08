###############################################################
# Comparison
###############################################################

(==)(a::MultiValue,b::MultiValue) = a.data == b.data
(≈)(a::MultiValue,b::MultiValue) = isapprox(get_array(a), get_array(b))
(≈)(a::MultiValue{S,T1,N,0} where T1,b::MultiValue{S,T2,N,0} where T2) where {S,N} = true

function (≈)(
  a::AbstractArray{<:MultiValue}, b::AbstractArray{<:MultiValue})
  if size(a) != size(b); return false; end
  for (ai,bi) in zip(a,b)
    if !(ai≈bi); return false; end
  end
  true
end

function isless(a::MultiValue{Tuple{L}},b::MultiValue{Tuple{L}}) where L
  for d in L:-1:1
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

    function ($op)(a::MultiValue{S},b::MultiValue{S})  where S
      r = broadcast(($op), a.data, b.data)
      T = change_eltype(a,eltype(r))
      T(r)
    end

    function ($op)(a::TensorValue,b::SymTensorValue)
      @notimplemented
    end

    function ($op)(a::SymTensorValue,b::TensorValue)
      @notimplemented
    end

  end
end

###############################################################
# Matrix Division
###############################################################

function (\)(a::MultiValue{Tuple{D,D}} where D, b::MultiValue)
  r = get_array(a) \ get_array(b)
  T = change_eltype(b,eltype(r))
  T(r)
end

###############################################################
# Operations with other numbers
###############################################################

for op in (:+,:-,:*)
  @eval begin
    function ($op)(a::MultiValue,b::Number)
        r = broadcast($op,a.data,b)
        T  = change_eltype(a,eltype(r))
        T(r)
    end

    function ($op)(a::Number,b::MultiValue)
        r = broadcast($op,a,b.data)
        T  = change_eltype(b,eltype(r))
        T(r)
    end
  end
end

function (/)(a::MultiValue,b::Number)
    r = broadcast(/,a.data,b)
    P  = change_eltype(a,eltype(r))
    P(r)
end

###############################################################
# Dot product (simple contraction)
###############################################################

function (*)(a::MultiValue, b::MultiValue)
  #msg = """
  #Method (*)(::$(typeof(a)),::$(typeof(b))) has been removed
  #Use simple contraction LinearAlgebra.⋅ (\cdot) or full contraction Gridap.⊙ (\odot) instead.
  #"""
  #error(msg)
  dot(a,b)
end

dot(a::MultiValue{Tuple{D}}, b::MultiValue{Tuple{D}}) where D = inner(a,b)

dot(a::MultiValue,b::MultiValue) = @notimplemented

@generated function dot(a::MultiValue{Tuple{D1}}, b::MultiValue{Tuple{D1,D2}}) where {D1,D2}
    ss = String[]
    for j in 1:D2
        s = join([ "a[$i]*b[$i,$j]+" for i in 1:D1])
        push!(ss,s[1:(end-1)]*", ")
    end
    str = join(ss)
    Meta.parse("VectorValue{$D2}($str)")
end

@generated function dot(a::MultiValue{Tuple{D1,D2}}, b::VectorValue{D1}) where {D1,D2}
    ss = String[]
    for j in 1:D2
        s = join([ "a[$j,$i]*b[$i]+" for i in 1:D1])
        push!(ss,s[1:(end-1)]*", ")
    end
    str = join(ss)
    Meta.parse("VectorValue{$D1}($str)")
end

@generated function dot(a::MultiValue{Tuple{D1,D3}}, b::MultiValue{Tuple{D3,D2}}) where {D1,D2,D3}
    ss = String[]
    for j in 1:D2
        for i in 1:D1
            s = join([ "a[$i,$k]*b[$k,$j]+" for k in 1:D3])
            push!(ss,s[1:(end-1)]*", ")
        end
    end
    str = join(ss)
    Meta.parse("TensorValue{$D1,$D2}(($str))")
end

###############################################################
# Inner product (full contraction)
###############################################################

inner(a::Real,b::Real) = a*b

function inner(a::MultiValue, b::MultiValue)
  @notimplemented
end

@generated function inner(a::MultiValue{S}, b::MultiValue{S}) where S
    str = join([" a[$i]*b[$i] +" for i in 1:length(a) ])
    Meta.parse(str[1:(end-1)])
end

@generated function inner(a::SymTensorValue{D}, b::SymTensorValue{D}) where D
  str = ""
  for i in 1:D
    for j in 1:D
      k = _2d_tensor_linear_index(D,i,j)
      str *= " a.data[$k]*b.data[$k] +"
    end
  end
  Meta.parse(str[1:(end-1)])
end

const ⊙ = inner

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

outer(a::MultiValue,b::MultiValue) = @notimplemented

@generated function outer(a::MultiValue{Tuple{D}},b::MultiValue{Tuple{Z}}) where {D,Z}
    str = join(["a[$i]*b[$j], " for j in 1:Z for i in 1:D])
    Meta.parse("TensorValue{$D,$Z}(($str))")
end

@generated function outer(a::MultiValue{Tuple{D}},b::MultiValue{Tuple{D1,D2}}) where {D,D1,D2}
  str = join(["a.array[$i]*b.array[$j,$k], "  for k in 1:D2 for j in 1:D1 for i in 1:D])
  Meta.parse("ThirdOrderTensorValue($str))")
end

const ⊗ = outer


###############################################################
# Linear Algebra
###############################################################

#TODO: write det and inv function for small specific cases.
det(a::MultiValue{Tuple{D1,D2}}) where {D1,D2} = det(get_array(a))
inv(a::MultiValue{Tuple{D1,D2}}) where {D1,D2} = TensorValue(inv(get_array(a)))

###############################################################
# Measure
###############################################################

"""
"""
meas(a::MultiValue{Tuple{D}}) where D = sqrt(inner(a,a))
meas(a::MultiValue{Tuple{D,D}}) where D = abs(det(a))

function meas(v::MultiValue{Tuple{1,2}})
  n1 = v[1,2]
  n2 = -1*v[1,1]
  n = VectorValue(n1,n2)
  sqrt(n*n)
end

function meas(v::MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1,n2,n3)
  sqrt(n*n)
end

@inline norm(u::MultiValue{Tuple{D}}) where D = sqrt(inner(u,u))
@inline norm(u::VectorValue{Tuple{0}}) = sqrt(zero(T))

###############################################################
# conj
###############################################################

function conj(a::T) where {T<:MultiValue}
  r = map(conj, a.data)
  T(r)
end

###############################################################
# Trace
###############################################################

@generated function tr(v::MultiValue{Tuple{D,D}}) where D
    str = join([" v[$i,$i] +" for i in 1:D ])
    Meta.parse(str[1:(end-1)])
end

@generated function tr(v::MultiValue{Tuple{A,A,B}}) where {A,B}
  lis = LinearIndices((A,A,B))
  str = ""
  for k in 1:B
    for i in 1:A
      if i !=1
        str *= " + "
      end
      p = lis[i,i,k]
      str *= " v.data[$p]"
    end
    str *= ", "
  end
  Meta.parse("VectorValue($str)")
end

###############################################################
# Adjoint and transpose
###############################################################

adjoint(a::MultiValue{Tuple{D,D}}) where D = @notimplemented
transpose(a::MultiValue{Tuple{D,D}}) where D = @notimplemented

@generated function adjoint(a::TensorValue{D1,D2}) where {D1,D2}
  str = ""
  for i in 1:D1
    for j in 1:D2
      k = (j-1)*D1 + i
      str *= "conj(a.data[$k]), "
    end
  end
  Meta.parse("TensorValue{D2,D1}($str)")
end

@generated function transpose(a::TensorValue{D1,D2}) where {D1,D2}
  str = ""
  for i in 1:D1
    for j in 1:D2
      k = (j-1)*D1 + i
      str *= "a.data[$k], "
    end
  end
  Meta.parse("TensorValue{D2,D1}($str)")
end

@inline function adjoint(a::TensorValue{D1,D2,T}) where {D1,D2,T<:Real}
  transpose(a)
end

adjoint(a::SymTensorValue) = conj(a)

@inline adjoint(a::SymTensorValue{D,T} where {D,T<:Real}) = transpose(a)

transpose(a::SymTensorValue) = a

###############################################################
# Symmetric part
###############################################################

"""
"""
@generated function symmetric_part(v::MultiValue{Tuple{D,D}}) where D
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

for op in (:symmetric_part,)
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

