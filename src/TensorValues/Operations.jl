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
      T = _eltype($op,r,a,b)
      M = change_eltype(a,T)
      M(r)
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
        T = _eltype($op,r,a,b)
        M  = change_eltype(a,T)
        M(r)
    end

    function ($op)(a::Number,b::MultiValue)
        r = broadcast($op,a,b.data)
        T = _eltype($op,r,a,b)
        M  = change_eltype(b,T)
        M(r)
    end
  end
end

function (/)(a::MultiValue,b::Number)
    r = broadcast(/,a.data,b)
    T = _eltype(/,r,a,b)
    P  = change_eltype(a,T)
    P(r)
end

@inline function _eltype(op,r,a,b)
  eltype(r)
end

@inline function _eltype(op,r::Tuple{},a,b)
  typeof(op(zero(eltype(a)),zero(eltype(b))))
end

###############################################################
# Dot product (simple contraction)
###############################################################

function (*)(a::MultiValue, b::MultiValue)
  msg = """
  Method (*)(::$(typeof(a)),::$(typeof(b))) has been removed.
  Depending the case, use simple contraction dot aka ⋅ (\\cdot) or full contraction inner aka ⊙ (\\odot) instead.
  """
  error(msg)
  #dot(a,b)
end

dot(a::MultiValue{Tuple{D}}, b::MultiValue{Tuple{D}}) where D = inner(a,b)

dot(a::MultiValue,b::MultiValue) = @notimplemented

@generated function dot(a::A,b::B) where {A<:MultiValue{Tuple{D1}},B<:MultiValue{Tuple{D1,D2}}} where {D1,D2}
    ss = String[]
    for j in 1:D2
      s = ""
      for i in 1:D1
        ak = data_index(A,i)
        bk = data_index(B,i,j)
        s *= "a.data[$ak]*b.data[$bk]+"
      end
        push!(ss,s[1:(end-1)]*", ")
    end
    str = join(ss)
    Meta.parse("VectorValue{$D2}($str)")
end

@generated function dot(a::A,b::B) where {A<:MultiValue{Tuple{D1,D2}},B<:MultiValue{Tuple{D2}}} where {D1,D2}
    ss = String[]
    for i in 1:D1
      s = ""
      for j in 1:D2
        ak = data_index(A,i,j)
        bk = data_index(B,j)
        s *= "a.data[$ak]*b.data[$bk]+"
      end
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

# a_ij = b_ijk*c_k
@generated function dot(a::A, b::B) where {A<:MultiValue{Tuple{D1,D2,D3}},B<:MultiValue{Tuple{D3}}} where {D1,D2,D3}
  ss = String[]
  for j in 1:D2
    for i in 1:D1
      s = join([ "a[$i,$j,$k]*b[$k]+" for k in 1:D3])
      push!(ss,s[1:(end-1)]*", ")
    end
  end
  str = join(ss)
  Meta.parse("TensorValue{$D1,$D2}($str)")
end

# a_ij = c_k*b_kij
@generated function dot(a::A, b::B) where {A<:MultiValue{Tuple{D1}},B<:MultiValue{Tuple{D1,D2,D3}}} where {D1,D2,D3}
  ss = String[]
  for k in 1:D3
    for j in 1:D2
      s = join([ "a[$i]*b[$i,$j,$k]+" for i in 1:D1])
      push!(ss,s[1:(end-1)]*", ")
    end
  end
  str = join(ss)
  Meta.parse("TensorValue{$D2,$D3}($str)")
end

# a_ilm = b_ij*c_jlm
@generated function dot(a::A,b::B) where {A<:MultiValue{Tuple{D,D}},B<:ThirdOrderTensorValue{D}} where D
  ss = String[]
  for m in 1:D
    for l in 1:D
      for i in 1:D
        s = join([ "a[$i,$j]*b[$j,$l,$m]+" for j in 1:D])
        push!(ss,s[1:(end-1)]*", ")
      end
    end
  end
  str = join(ss)
  Meta.parse("ThirdOrderTensorValue{$D}($str)")
end

const ⋅¹ = dot

###############################################################
# Inner product (full contraction)
###############################################################

inner(a::Number,b::Number) = a*b

function inner(a::Complex,b::Complex)
  @notimplemented "What do we want here a*b or a'*b ?"
end

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
      k = data_index(a,i,j)
      str *= " a.data[$k]*b.data[$k] +"
    end
  end
  Meta.parse(str[1:(end-1)])
end

@generated function inner(a::SymFourthOrderTensorValue{D}, b::SymTensorValue{D}) where D
  str = ""
  for i in 1:D
    for j in i:D
      s = ""
      for k in 1:D
        for l in 1:D
          ak = data_index(a,i,j,k,l)
          bk = data_index(b,k,l)
          s *= " a.data[$ak]*b.data[$bk] +"
        end
      end
      str *= s[1:(end-1)]*", "
    end
  end
  Meta.parse("SymTensorValue{D}($str)")
end

function inner(a::SymFourthOrderTensorValue{D},b::MultiValue{Tuple{D,D}}) where D
  inner(a,symmetric_part(b))
end

const ⊙ = inner

###############################################################
# Double Contractions w/ products
###############################################################

# a_i = b_ijk*c_jk
@generated function double_contraction(a::A, b::B) where {A<:MultiValue{Tuple{D1,D2,D3}},B<:MultiValue{Tuple{D2,D3}}} where {D1,D2,D3}
  ss = String[]
  for i in 1:D1
    s = join([ "a[$i,$j,$k]*b[$j,$k]+" for j in 1:D2 for k in 1:D3])
    push!(ss,s[1:(end-1)]*", ")
  end
  str = join(ss)
  Meta.parse("VectorValue{$D1}(($str))")
end

# a_ijpm = b_ijkl*c_klpm (3D)
@generated function double_contraction(a::A, b::B) where {A<:SymFourthOrderTensorValue{3},B<:SymFourthOrderTensorValue{3}}

  Sym4TensorIndexing = [1111, 1121, 1131, 1122, 1132, 1133, 2111, 2121, 2131, 2122, 2132, 2133,
                        3111, 3121, 3131, 3122, 3132, 3133, 2211, 2221, 2231, 2222, 2232, 2233,
                        2311, 2321, 2331, 2322, 2332, 2333, 3311, 3321, 3331, 3322, 3332, 3333]
  ss = String[]
  for off_index in Sym4TensorIndexing
    i = parse(Int,string(off_index)[1]); j = parse(Int,string(off_index)[2]);
    m = parse(Int,string(off_index)[3]); p = parse(Int,string(off_index)[4]);
    s = join([ "a[$i,$j,$k,$l]*b[$k,$l,$m,$p]+" for k in 1:3 for l in 1:3])
    push!(ss,s[1:(end-1)]*", ")
  end
  str = join(ss)
  Meta.parse("SymFourthOrderTensorValue{3}($str)")
end

# a_ijpm = b_ijkl*c_klpm (general case)
@generated function double_contraction(a::SymFourthOrderTensorValue{D}, b::SymFourthOrderTensorValue{D}) where D
  str = ""
  for j in 1:D
    for i in j:D
      for m in 1:D
        for p in m:D
          s = ""
          for k in 1:D
            for l in 1:D
              s *= " a[$i,$j,$k,$l]*b[$k,$l,$p,$m] +"
            end
          end
          str *= s[1:(end-1)]*", "
        end
      end
    end
  end
  Meta.parse("SymFourthOrderTensorValue{D}($str)")
end

# a_ilm = b_ijk*c_jklm
@generated function double_contraction(a::A,b::B) where {A<:ThirdOrderTensorValue{D},B<:SymFourthOrderTensorValue{D}} where D
  ss = String[]
  for m in 1:D
    for l in 1:D
      for i in 1:D
        s = join([ "a[$i,$j,$k]*b[$j,$k,$l,$m]+" for j in 1:D for k in 1:D])
        push!(ss,s[1:(end-1)]*", ")
      end
    end
  end
  str = join(ss)
  Meta.parse("ThirdOrderTensorValue{$D}($str)")
end

# a_il = b_ijk*c_jkl
@generated function double_contraction(a::A,b::B) where {A<:ThirdOrderTensorValue{D},B<:ThirdOrderTensorValue{D}} where D
  ss = String[]
  for l in 1:D
    for i in 1:D
      s = join([ "a[$i,$j,$k]*b[$j,$k,$l]+" for j in 1:D for k in 1:D])
      push!(ss,s[1:(end-1)]*", ")
    end
  end
  str = join(ss)
  Meta.parse("TensorValue{$D}($str)")
end

const ⋅² = double_contraction

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
outer(a::Number,b::Number) = a*b

outer(a::MultiValue,b::Number) = a*b
outer(a::Number,b::MultiValue) = a*b

function outer(a::MultiValue,b::MultiValue)
   @notimplemented
end

@generated function outer(a::MultiValue{Tuple{D}},b::MultiValue{Tuple{Z}}) where {D,Z}
    str = join(["a[$i]*b[$j], " for j in 1:Z for i in 1:D])
    Meta.parse("TensorValue{$D,$Z}($str)")
end

function outer(a::VectorValue{0,Ta},b::VectorValue{1,Tb}) where {Ta,Tb}
  T = promote_type(Ta,Tb)
  TensorValue{0,1,T}()
end

function outer(a::VectorValue{0,Ta},b::Tb) where {Ta,Tb<:Number}
  T = promote_type(Ta,Tb)
  VectorValue{0,T}()
end

@generated function outer(a::MultiValue{Tuple{D}},b::MultiValue{Tuple{D1,D2}}) where {D,D1,D2}
  str = join(["a[$i]*b[$j,$k], "  for k in 1:D2 for j in 1:D1 for i in 1:D])
  Meta.parse("ThirdOrderTensorValue{D,D1,D2}($str)")
end

@generated function outer(a::SymTensorValue{D},b::SymTensorValue{D}) where D
  str = ""
  for i in 1:D
    for j in i:D
      ak = data_index(a,i,j)
      for k in 1:D
        for l in k:D
          bk = data_index(b,k,l)
          str *= "a.data[$ak]*b.data[$bk], "
        end
      end
    end
  end
  Meta.parse("SymFourthOrderTensorValue{D}($str)")
end

const ⊗ = outer

###############################################################
# Cross Product
###############################################################

function cross(a::MultiValue{Tuple{3}}, b::MultiValue{Tuple{3}})
  VectorValue{3}(a[2]b[3]-a[3]b[2], a[3]b[1]-a[1]b[3], a[1]b[2]-a[2]b[1])
end

function cross(a::MultiValue{Tuple{2}}, b::MultiValue{Tuple{2}})
  a[1]b[2]-a[2]b[1]
end

cross(a::MultiValue,b::MultiValue) = error("Cross product only defined for R2 and R3 vectors")

###############################################################
# Linear Algebra
###############################################################

det(a::MultiValue{Tuple{D1,D2}}) where {D1,D2} = det(get_array(a))

det(a::MultiValue{Tuple{1,1}}) = a[1]

function det(a::MultiValue{Tuple{2,2}})
  a_11 = a[1,1]; a_12 = a[1,2]
  a_21 = a[2,1]; a_22 = a[2,2]
  a_11*a_22 - a_12*a_21 
end

function det(a::MultiValue{Tuple{3,3}})
  a_11 = a[1,1]; a_12 = a[1,2]; a_13 = a[1,3]
  a_21 = a[2,1]; a_22 = a[2,2]; a_23 = a[2,3]
  a_31 = a[3,1]; a_32 = a[3,2]; a_33 = a[3,3]
  a_11*a_22*a_33 + a_12*a_23*a_31 + a_13*a_21*a_32 -
    (a_11*a_23*a_32 + a_12*a_21*a_33 + a_13*a_22*a_31)
end

inv(a::MultiValue{Tuple{D1,D2}}) where {D1,D2} = TensorValue(inv(get_array(a)))

function inv(a::MultiValue{Tuple{1,1}})
  r = 1/a[1]
  T = change_eltype(a,typeof(r))
  T(r)
end

function inv(a::MultiValue{Tuple{2,2}})
 c = 1/det(a)
 data = (a[2,2]*c, -a[2,1]*c, -a[1,2]*c,  a[1,1]*c)
 TensorValue{2}(data)
end

function inv(a::MultiValue{Tuple{3,3}})
 a_11 = a[1,1]; a_12 = a[1,2]; a_13 = a[1,3]
 a_21 = a[2,1]; a_22 = a[2,2]; a_23 = a[2,3]
 a_31 = a[3,1]; a_32 = a[3,2]; a_33 = a[3,3]
 c = 1/det(a)
 data = (
    ( a_22*a_33 - a_23*a_32 )*c,
   -( a_21*a_33 - a_23*a_31 )*c,
    ( a_21*a_32 - a_22*a_31 )*c,
   -( a_12*a_33 - a_13*a_32 )*c,
    ( a_11*a_33 - a_13*a_31 )*c,
   -( a_11*a_32 - a_12*a_31 )*c,
    ( a_12*a_23 - a_13*a_22 )*c,
   -( a_11*a_23 - a_13*a_21 )*c,
    ( a_11*a_22 - a_12*a_21 )*c)
 TensorValue{3}(data)
end

###############################################################
# Measure
###############################################################

"""
"""
meas(a::MultiValue{Tuple{D}}) where D = sqrt(inner(a,a))
meas(a::MultiValue{Tuple{D,D}}) where D = abs(det(a))
meas(a::TensorValue{0,1,T}) where T = one(T)

function meas(v::MultiValue{Tuple{1,2}})
  n1 = v[1,2]
  n2 = -1*v[1,1]
  n = VectorValue(n1,n2)
  sqrt(n ⋅ n)
end

function meas(v::MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1,n2,n3)
  sqrt(n ⋅ n)
end

@inline norm(u::MultiValue{Tuple{D}}) where D = sqrt(inner(u,u))
@inline norm(u::MultiValue{Tuple{0},T}) where T = sqrt(zero(T))

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
        for i in j:D
            str *= "0.5*v[$i,$j] + 0.5*v[$j,$i], "
        end
    end
    str *= ")"
    Meta.parse("SymTensorValue{D}($str)")
end

###############################################################
# Define new operations for Gridap types
###############################################################

#for op in (:symmetric_part,)
#    @eval begin
#        ($op)(a::GridapType) = operate($op,a)
#    end
#end
#
#for op in (:inner,:outer,:double_contraction)#,:(:))
#    @eval begin
#        ($op)(a::GridapType,b::GridapType) = operate($op,a,b)
#        ($op)(a::GridapType,b::Number)     = operate($op,a,b)
#        ($op)(a::Number,    b::GridapType) = operate($op,a,b)
#        ($op)(a::GridapType,b::Function)   = operate($op,a,b)
#        ($op)(a::Function,  b::GridapType) = operate($op,a,b)
#    end
#end



