###############################################################
# Comparison
###############################################################

(==)(a::MultiValue,b::MultiValue) = false
(==)(a::MultiValue{S},b::MultiValue{S}) where {S} = a.data == b.data
(≈)(a::MultiValue{S},b::MultiValue{S}) where {S} = isapprox(get_array(a), get_array(b))
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
    if isless(a[d], b[d])
      return true
    elseif isless(b[d], a[d])
      return false
    else
      continue
    end
  end
  false
end

isless(a::Number,b::MultiValue) = all(isless.(a, b.data))
isless(a::MultiValue,b::MultiValue) = @unreachable "Comparison is not defined between tensor of order greater than 1"

###############################################################
# Addition / subtraction
###############################################################

Base.iszero(a::MultiValue) = all(iszero.(a.data))

for op in (:+,:-)
  @eval begin

    function ($op)(a::T) where {T<:MultiValue}
      r = map($op, a.data)
      T(r)
    end

    function ($op)(a::MultiValue,b::MultiValue)
      @notimplemented "Not implemented or undefined operation \"$($op)\" on MultiValues of these shapes"
    end

    function ($op)(a::MultiValue{S},b::MultiValue{S})  where S
      r = map(($op), a.data, b.data)
      T = _eltype($op,r,a,b)
      M = change_eltype(a,T)
      M(r)
    end

    function ($op)(a::TensorValue{D,D},b::SymTensorValue{D}) where D
      map(($op), a, TensorValue(get_array(b)))
    end

    function ($op)(a::SymTensorValue{D},b::TensorValue{D,D}) where D
      map(($op), TensorValue(get_array(a)), b)
    end

    function ($op)(a::TensorValue{D,D},b::SymTracelessTensorValue{D}) where D
      map(($op), a, TensorValue(get_array(b)))
    end

    function ($op)(a::SymTracelessTensorValue{D},b::TensorValue{D,D}) where D
      map(($op), TensorValue(get_array(a)), b)
    end

    function ($op)(a::SymTracelessTensorValue{D},b::SymTensorValue{D}) where D
      r = map(($op), a.data, b.data)
      T = _eltype($op,r,a,b)
      M = change_eltype(b,T)
      M(r)
    end

    function ($op)(a::SymTensorValue{D},b::SymTracelessTensorValue{D}) where D
      r = map(($op), a.data, b.data)
      T = _eltype($op,r,a,b)
      M = change_eltype(a,T)
      M(r)
    end

    function ($op)(a::SymTracelessTensorValue)
      r = map($op, a.data[1:end-1])
      typeof(a)(r)
    end

    function ($op)(a::SymTracelessTensorValue{D},b::SymTracelessTensorValue{D}) where D
      r = map(($op), a.data[1:end-1], b.data[1:end-1])
      T = _eltype($op,r,a,b)
      M = change_eltype(a,T)
      M(r)
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

@generated function _bc(f,a::NTuple{N},b::Number) where N
  s = "("
  for i in 1:N
    s *= "f(a[$i],b), "
  end
  s *= ")"
  Meta.parse(s)
end

@generated function _bc(f,b::Number,a::NTuple{N}) where N
  s = "("
  for i in 1:N
    s *= "f(b,a[$i]), "
  end
  s *= ")"
  Meta.parse(s)
end

for op in (:+,:-,:*)
  @eval begin
    function ($op)(a::MultiValue,b::Number)
      r = _bc($op,a.data,b)
      T = _eltype($op,r,a,b)
      M  = change_eltype(a,T)
      M(r)
    end

    function ($op)(a::Number,b::MultiValue)
      r = _bc($op,a,b.data)
      T = _eltype($op,r,a,b)
      M  = change_eltype(b,T)
      M(r)
    end
  end
end

function (*)(a::Number,b::SymTracelessTensorValue)
  r = _bc(*,a,b.data[1:end-1])
  T = _eltype(*,r,a,b)
  M  = change_eltype(b,T)
  M(r)
end

function (*)(a::SymTracelessTensorValue,b::Number)
  b*a
end

function (/)(a::MultiValue,b::Number)
  r = _bc(/,a.data,b)
  T = _eltype(/,r,a,b)
  P  = change_eltype(a,T)
  P(r)
end

function (/)(a::SymTracelessTensorValue,b::Number)
  r = _bc(/,a.data[1:end-1],b)
  T = _eltype(/,r,a,b)
  M  = change_eltype(a,T)
  M(r)
end

const _err =  " with number is undefined for traceless tensors"
function +(::SymTracelessTensorValue,::Number)     error("Addition"   *_err) end
function -(::SymTracelessTensorValue,::Number)     error("Subtraction"*_err) end
function +(::Number,::SymTracelessTensorValue)     error("Addition"   *_err) end
function -(::Number,::SymTracelessTensorValue)     error("Subtraction"*_err) end
function +(::SymTracelessTensorValue,::MultiValue) error("Addition"   *_err) end
function -(::SymTracelessTensorValue,::MultiValue) error("Subtraction"*_err) end
function +(::MultiValue,::SymTracelessTensorValue) error("Addition"   *_err) end
function -(::MultiValue,::SymTracelessTensorValue) error("Subtraction"*_err) end

@inline function _eltype(op,r,a...)
  eltype(r)
end

@inline function _eltype(op,r::Tuple{},a...)
  typeof(reduce(op,zero.(eltype.(a))))
end

@inline function _eltype(op,r,a,b)
  eltype(r)
end

@inline function _eltype(op,r::Tuple{},a,b)
  typeof(op(zero(eltype(a)),zero(eltype(b))))
end

@inline function _eltype(op,r,a)
  eltype(r)
end

@inline function _eltype(op,r::Tuple{},a)
  typeof(op(zero(eltype(a))))
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
end

# Resolution of silly method ambiguity
const _msg =  "Use use simple contraction dot aka ⋅ (\\cdot) or full contraction inner aka ⊙ (\\odot)"
function *(::MultiValue,::SymTracelessTensorValue) @unreachable _msg end
function *(::SymTracelessTensorValue,::MultiValue) @unreachable _msg end
function *(::SymTracelessTensorValue,::AbstractSymTensorValue) @unreachable _msg end
function *(::SymTracelessTensorValue,::SymTracelessTensorValue) @unreachable _msg end

dot(a::MultiValue{Tuple{D}}, b::MultiValue{Tuple{D}}) where D = inner(a,b)

"""
    dot(a::MultiValue{Tuple{...,D}}, b::MultiValue{Tuple{D,...}})
    a ⋅¹ b
    a ⋅ b

Inner product of two tensors `a` and `b`, that is the single contraction of the last index of `a` with the first index of `b`. The corresponding dimensions `D` must match. No symmetry is preserved.
"""
dot(a::MultiValue,b::MultiValue) = @notimplemented

@generated function dot(a::MultiValue{Tuple{D1},Ta},b::MultiValue{Tuple{D1,D2},Tb}) where {D1,D2,Ta,Tb}
  iszero(length(b)) && return :( zero(VectorValue{D2,$(promote_type(Ta,Tb))}) )
  ss = String[]
  for j in 1:D2
    s = ""
    for i in 1:D1
      s *= "a[$i]*b[$i,$j]+"
    end
    push!(ss,s[1:(end-1)]*", ")
  end
  str = join(ss)
  Meta.parse("VectorValue{$D2}($str)")
end

@generated function dot(a::MultiValue{Tuple{D1,D2},Ta},b::MultiValue{Tuple{D2},Tb}) where {D1,D2,Ta,Tb}
  iszero(length(a)) && return :( zero(VectorValue{D1,$(promote_type(Ta,Tb))}) )
  ss = String[]
  for i in 1:D1
    s = ""
    for j in 1:D2
      s *= "a[$i,$j]*b[$j]+"
    end
      push!(ss,s[1:(end-1)]*", ")
  end
  str = join(ss)
  Meta.parse("VectorValue{$D1}($str)")
end

@generated function dot(a::MultiValue{Tuple{D1,D3},Ta}, b::MultiValue{Tuple{D3,D2},Tb}) where {D1,D2,D3,Ta,Tb}
  (iszero(length(a)) || iszero(length(b))) && return :( zero(TensorValue{D1,D2,$(promote_type(Ta,Tb))}) )
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
@generated function dot(a::MultiValue{Tuple{D1,D2,D3},Ta}, b::MultiValue{Tuple{D3},Tb}) where {D1,D2,D3,Ta,Tb}
  iszero(length(a)) && return :( zero(TensorValue{D1,D2,$(promote_type(Ta,Tb))}) )
  T = promote_type(Ta,Tb)
  (iszero(D1) || iszero(D2)) && return :( TensorValue{D1,D2,$T}() )
  iszero(D3) && return :( zero(TensorValue{D1,D2,$T}) )
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

# a_ijl = b_ijk*c_kl
@generated function dot(a::MultiValue{Tuple{D1,D2,D3},Ta}, b::MultiValue{Tuple{D3,D4},Tb}) where {D1,D2,D3,D4,Ta,Tb}
  (iszero(length(a)) || iszero(length(b))) && return :(
    zero(ThirdOrderTensorValue{D1,D2,D4,$(promote_type(Ta,Tb))})
  )
  ss = String[]
  for l in 1:D4
    for j in 1:D2
      for i in 1:D1
        s = join([ "a[$i,$j,$k]*b[$k,$l]+" for k in 1:D3])
        push!(ss,s[1:(end-1)]*", ")
      end
    end
  end
  str = join(ss)
  Meta.parse("ThirdOrderTensorValue{$D1,$D2,$D4}($str)")
end

# a_ij = c_k*b_kij
@generated function dot(a::MultiValue{Tuple{D1},Ta}, b::MultiValue{Tuple{D1,D2,D3},Tb}) where {D1,D2,D3,Ta,Tb}
  iszero(length(b)) && return :( zero(TensorValue{D2,D3,$(promote_type(Ta,Tb))}) )
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
@generated function dot(a::MultiValue{Tuple{D1,D2},Ta},b::MultiValue{Tuple{D2,D3,D4},Tb}) where  {D1,D2,D3,D4,Ta,Tb}
  (iszero(length(a)) || iszero(length(b))) && return :(
    zero(ThirdOrderTensorValue{D1,D3,D4,$(promote_type(Ta,Tb))})
  )
  ss = String[]
  for m in 1:D4
    for l in 1:D3
      for i in 1:D1
        s = join([ "a[$i,$j]*b[$j,$l,$m]+" for j in 1:D2])
        push!(ss,s[1:(end-1)]*", ")
      end
    end
  end
  str = join(ss)
  Meta.parse("ThirdOrderTensorValue{$D1,$D3,$D4}($str)")
end

const ⋅¹ = dot

###############################################################
# Inner product (full contraction)
###############################################################

inner(a::Number,b::Number) = a*b

"""
    inner(a::MultiValue{S}, b::MultiValue{S}) -> scalar
    a ⊙ b

Inner product of two tensors, that is the full contraction along each indices. The size `S` of `a` and `b` must match.
"""
function inner(a::MultiValue, b::MultiValue)
  @notimplemented "Sizes of tensors must match."
end

@generated function inner(a::MultiValue{S,Ta}, b::MultiValue{S,Tb}) where {S,Ta,Tb}
  iszero(length(a)) && return :( zero($(promote_type(Ta,Tb))) )
  str = join([" a[$i]*b[$i] +" for i in 1:length(a) ])
  Meta.parse(str[1:(end-1)])
end

@generated function inner(a::AbstractSymTensorValue{D,Ta}, b::AbstractSymTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :( zero($(promote_type(Ta,Tb))) )
  str = ""
  for i in 1:D
    str *= "+ a[$i,$i]*b[$i,$i]"
  end
  str *= " + 2*("
  for i in 1:D
    for j in i+1:D
      str *= "+ a[$i,$j]*b[$i,$j]"
    end
  end
  str *= ")"
  Meta.parse(str)
end

function inner(a::MultiValue{Tuple{D,D,D,D}}, b::MultiValue{Tuple{D,D,D,D}}) where D
  double_contraction(a,b)
end

function inner(a::MultiValue{Tuple{D,D,D,D}}, b::MultiValue{Tuple{D,D}}) where D
  double_contraction(a,b)
end

function inner(a::MultiValue{Tuple{D,D}}, b::MultiValue{Tuple{D,D,D,D}}) where D
  double_contraction(a,b)
end

const ⊙ = inner

###############################################################
# Double Contractions w/ products
###############################################################

"""
    double_contraction(a::MultiValue{Tuple{...,D,E}}, b::MultiValue{Tuple{D,E,...})
    a ⋅² b

Double contraction of two tensors `a` and `b`, along the two last indices of `a`
and two first of `b`. The corresponding dimensions `D` and `E` must match, the
contraction order is chosen to be consistent with the inner product of second
order tensors.

The `double_contraction` between second- and/or fourth-order symmetric tensors
preserves the symmetry (returns a symmetric tensor type).
"""
function double_contraction(a::MultiValue{S1}, b::MultiValue{S2}) where {S1<:Tuple,S2<:Tuple}
  L1, L2 = length(S1.types), length(S2.types)
  if L1<2 || L2<2
    @unreachable "Double contraction is only define for tensors of order more than 2, got $L1 and $L2."
  end

  D1, E1, D2, E2 = S1.types[end-1], S1.types[end],  S2.types[1], S2.types[2]
  if D1 != D2 || E1 != E2
    throw(DimensionMismatch("the last two dimensions of the first argument must match the first two of the second argument, got ($D1,$E1) ≠ ($D2,$E2)."))
  end
  @notimplemented
end

# c_i = a_ij*b_ij
function double_contraction(a::MultiValue{S}, b::MultiValue{S}) where {S<:Tuple{D1,D2}} where {D1,D2}
  inner(a,b)
end

# c_i = a_ijk*b_jk
@generated function double_contraction(a::MultiValue{Tuple{D1,D2,D3},Ta}, b::MultiValue{Tuple{D2,D3},Tb})  where {D1,D2,D3,Ta,Tb}
  iszero(length(a)) && return :( zero(VectorValue{D1,$(promote_type(Ta,Tb))}) )
  ss = String[]
  for i in 1:D1
    s = join([ "a[$i,$j,$k]*b[$j,$k]+" for j in 1:D2 for k in 1:D3])
    push!(ss,s[1:(end-1)]*", ")
  end
  str = join(ss)
  Meta.parse("VectorValue{$D1}(($str))")
end

# c_k = a_ij*b_ijk
@generated function double_contraction(a::MultiValue{Tuple{D1,D2},Ta}, b::MultiValue{Tuple{D1,D2,D3},Tb})  where {D1,D2,D3,Ta,Tb}
  iszero(length(b)) && return :( zero(VectorValue{D3,$(promote_type(Ta,Tb))}) )
  ss = String[]
  for k in 1:D3
    s = join([ "a[$i,$j]*b[$i,$j,$k]+" for i in 1:D1 for j in 1:D2])
    push!(ss,s[1:(end-1)]*", ")
  end
  str = join(ss)
  Meta.parse("VectorValue{$D3}(($str))")
end

# c_ijpm = a_ijkl*b_klpm (3D)
@generated function double_contraction(a::SymFourthOrderTensorValue{3}, b::SymFourthOrderTensorValue{3})

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

# c_ijpm = a_ijkl*b_klpm (general case)
@generated function double_contraction(a::SymFourthOrderTensorValue{D,Ta}, b::SymFourthOrderTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :( SymFourthOrderTensorValue{0,$(promote_type(Ta,Tb))}() )
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

# c_ilm = a_ijk*b_jklm
@generated function double_contraction(a::ThirdOrderTensorValue{D1,D,D,Ta},b::SymFourthOrderTensorValue{D,Tb}) where {D1,D,Ta,Tb}
  iszero(length(a)) && return :( zero(ThirdOrderTensorValue{D1,D,D,$(promote_type(Ta,Tb))}) )
  ss = String[]
  for m in 1:D
    for l in 1:D
      for i in 1:D1
        s = join([ "a[$i,$j,$k]*b[$j,$k,$l,$m]+" for j in 1:D for k in 1:D])
        push!(ss,s[1:(end-1)]*", ")
      end
    end
  end
  str = join(ss)
  Meta.parse("ThirdOrderTensorValue{$D1,$D,$D}($str)")
end

# c_ij = a_ijkl*b_kl
@generated function double_contraction(a::SymFourthOrderTensorValue{D,Ta}, b::AbstractSymTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :( zero(SymTensorValue{D,$(promote_type(Ta,Tb))}) )
  str = ""
  for i in 1:D
    for j in i:D
      for k in 1:D
        str *= "+ a[$i,$j,$k,$k]*b[$k,$k]"
      end
      str *= " + 2*("
      for k in 1:D
        for l in k+1:D
          str *= "+ a[$i,$j,$k,$l]*b[$k,$l]"
        end
      end
      str *= "), "
    end
  end
  Meta.parse("SymTensorValue{D}($str)")
end

# c_kl = a_ij*b_ijkl
@generated function double_contraction(a::AbstractSymTensorValue{D,Ta}, b::SymFourthOrderTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :( zero(SymTensorValue{D,$(promote_type(Ta,Tb))}) )
  str = ""
  for k in 1:D
    for l in k:D
      for i in 1:D
        str *= "+ a[$i,$i]*b[$i,$i,$k,$l]"
      end
      str *= " + 2*("
      for i in 1:D
        for j in i+1:D
          str *= "+ a[$i,$j]*b[$i,$j,$k,$l]"
        end
      end
      str *= "), "
    end
  end
  Meta.parse("SymTensorValue{D}($str)")
end

# c_ij = a_ijkl*b_kl
function double_contraction(a::SymFourthOrderTensorValue{D}, b::MultiValue{Tuple{D,D}}) where D
  double_contraction(a,symmetric_part(b))
end

# c_kl = a_ij*b_ijkl
function double_contraction(a::MultiValue{Tuple{D,D}}, b::SymFourthOrderTensorValue{D}) where D
  double_contraction(symmetric_part(a),b)
end


# c_il = a_ijk*b_jkl
@generated function double_contraction(a::ThirdOrderTensorValue{D1,D,E,Ta},b::ThirdOrderTensorValue{D,E,D2,Tb}) where {D1,D,E,D2,Ta,Tb}
  (iszero(length(a)) || iszero(length(b))) && return :(
    zero(TensorValue{D1,D2,$(promote_type(Ta,Tb))})
  )
  ss = String[]
  for l in 1:D2
    for i in 1:D1
      s = join([ "a[$i,$j,$k]*b[$j,$k,$l]+" for j in 1:D for k in 1:E])
      push!(ss,s[1:(end-1)]*", ")
    end
  end
  str = join(ss)
  Meta.parse("TensorValue{$D1,$D2}($str)")
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

outer(a::Number,b::Number) = a*b

outer(a::MultiValue,b::Number) = a*b
outer(a::Number,b::MultiValue) = a*b

"""
    outer(a,b)
    a ⊗ b

Outer product (or tensor-product) of two `Number`s and/or `MultiValue`s, that is
`(a⊗b)[i₁,...,iₙ,j₁,...,jₙ] = a[i₁,...,iₙ]*b[j₁,...,jₙ]`. This falls back to standard
multiplication if `a` or `b` is a scalar.
"""
function outer(a::MultiValue,b::MultiValue)
   @notimplemented
end

@generated function outer(a::MultiValue{Tuple{D},Ta},b::MultiValue{Tuple{Z},Tb}) where {D,Z,Ta,Tb}
  (iszero(D) || iszero(Z)) && return :(
    zero(TensorValue{D,Z,$(promote_type(Ta,Tb))})
  )
  str = join(["a[$i]*b[$j], " for j in 1:Z for i in 1:D])
  Meta.parse("TensorValue{$D,$Z}($str)")
end

@generated function outer(a::MultiValue{Tuple{D},Ta},b::MultiValue{Tuple{D1,D2},Tb}) where {D,D1,D2,Ta,Tb}
  (iszero(D) || iszero(length(b))) && return :(
    zero(ThirdOrderTensorValue{D,D1,D2,$(promote_type(Ta,Tb))})
  )
  str = join(["a[$i]*b[$j,$k], "  for k in 1:D2 for j in 1:D1 for i in 1:D])
  Meta.parse("ThirdOrderTensorValue{D,D1,D2}($str)")
end

@generated function outer(a::MultiValue{Tuple{D1,D2},Ta},b::MultiValue{Tuple{D},Tb}) where {D,D1,D2,Ta,Tb}
  (iszero(length(a)) || iszero(D)) && return :(
    zero(ThirdOrderTensorValue{D1,D2,D,$(promote_type(Ta,Tb))})
  )
  str = join(["a[$i,$j]*b[$k], "  for k in 1:D for j in 1:D2 for i in 1:D1])
  Meta.parse("ThirdOrderTensorValue{D1,D2,D}($str)")
end

@generated function outer(a::AbstractSymTensorValue{D,Ta},b::AbstractSymTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :( zero(SymFourthOrderTensorValue{D,$(promote_type(Ta,Tb))}) )
  str = ""
  for i in 1:D
    for j in i:D
      for k in 1:D
        for l in k:D
          str *= "a[$i,$j]*b[$k,$l], "
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

"""
    cross(a::VectorValue{3}, b::VectorValue{3}) -> VectorValue{3}
    cross(a::VectorValue{2}, b::VectorValue{2}) -> Scalar
    a × b

Cross product of 2D and 3D vector.
"""
cross(a::MultiValue,b::MultiValue) = error("Cross product only defined for R2 and R3 vectors of same dimension")

###############################################################
# Linear Algebra
###############################################################

"""
    det(a::MultiValue{Tuple{D,D},T})

Determinent of square second order tensors.
"""
det(a::MultiValue{Tuple{D,D}}) where {D} = det(get_array(a))
det(a::MultiValue)= @unreachable "det undefined for this tensor shape: $(size(a))"

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

"""
    inv(a::MultiValue{Tuple{D,D}})

Inverse of a second order tensor.
"""
inv(a::MultiValue{Tuple{D,D}}) where D = TensorValue(inv(get_array(a)))

# this has better perf than the D=2,3 specialization below
inv(a::SymTracelessTensorValue{2}) = SymTracelessTensorValue(inv(get_array(a)))

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
    meas(a::MultiValue{Tuple{D}})
    meas(a::MultiValue{Tuple{1,D2}})

Euclidean norm of a vector.
"""
meas(a::MultiValue{Tuple{D}}) where D = sqrt(inner(a,a))

"""
    meas(J::MultiValue{Tuple{D1,D2}})

Returns the absolute `D1`-dimensional volume of the parallelepiped
formed by the rows of `J`, that is `sqrt(det(J⋅Jᵀ))`, or `abs(det(J))` if `D1`=`D2`.
This is used to compute the contribution of the Jacobian matrix `J` of a changes of variables in integrals.
"""
meas(a::MultiValue{Tuple{D,D}}) where D = abs(det(a))
#meas( ::TensorValue{0,D,T}) where {T,D} = one(T)
#meas( ::MultiValue{Tuple{0,0},T}) where {T} = one(T)

function meas(v::MultiValue{Tuple{1,D}}) where D
  t = VectorValue(v.data)
  meas(t)
end

function meas(v::MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1,n2,n3)
  meas(n)
end

function meas(Jt::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  J = transpose(Jt)
  sqrt(det(Jt⋅J))
end

"""
    norm(u::MultiValue{Tuple{D}})
    norm(u::MultiValue{Tuple{D1,D2}})

Euclidean (2-)norm of `u`, namely `sqrt(inner(u,u))`.
"""
@inline norm(u::MultiValue{Tuple{D},<:Real}) where D = sqrt(inner(u,u))
@inline norm(u::MultiValue{Tuple{D}}) where D = sqrt(real(inner(u,conj(u))))
@inline norm(u::MultiValue{Tuple{D1,D2},<:Real}) where {D1,D2} = sqrt(inner(u,u))
@inline norm(u::MultiValue{Tuple{D1,D2}}) where {D1,D2} = sqrt(real(inner(u,conj(u))))
@inline norm(u::MultiValue{Tuple{0},T}) where T = sqrt(zero(T))

###############################################################
# conj, real, imag
###############################################################

for op in (:conj,:real,:imag)
  @eval begin
    function ($op)(a::T) where {T<:MultiValue}
      r = map($op, a.data)
      T2 = _eltype($op,r,a)
      M  = change_eltype(a,T2)
      M(r)
    end

    function ($op)(a::T) where {T<:SymTracelessTensorValue}
      r = map($op, a.data)
      T2 = _eltype($op,r,a)
      M  = change_eltype(a,T2)
      M(r[1:end-1])
    end
  end
end

###############################################################
# Trace
###############################################################

"""
    tr(v::MultiValue{Tuple{D,D}})

Return the trace of a second order square tensor, defined by `Σᵢ vᵢᵢ` or 0 if `D`=0.
"""
@generated function tr(v::MultiValue{Tuple{D,D},T}) where {D,T}
  iszero(D) && return :(zero(T))
  str = join([" v[$i,$i] +" for i in 1:D ])
  Meta.parse(str[1:(end-1)])
end
tr(::SymTracelessTensorValue{D,T}) where {D,T} = zero(T)
tr(::MultiValue{Tuple{A,B}}) where {A,B} = throw(ArgumentError("Second order tensor is not square"))

"""
    tr(v::MultiValue{Tuple{D,D,D2}}) -> ::VectorValue{D2}

Return a vector of length `D2` of traces computed on the first two indices: `resⱼ = Σᵢ vᵢᵢⱼ`.
"""
@generated function tr(v::MultiValue{Tuple{A,A,B},T}) where {A,B,T}
  iszero(length(v)) && return :( zero(VectorValue{B,T}) )
  str = ""
  for k in 1:B
    for i in 1:A
      if i !=1
        str *= " + "
      end
      str *= " v[$i,$i,$k]"
    end
    str *= ", "
  end
  Meta.parse("VectorValue($str)")
end
tr(::MultiValue{Tuple{A,B,C}}) where {A,B,C} = throw(ArgumentError("First two dimensions are not iddentical"))

###############################################################
# Adjoint and transpose
###############################################################

adjoint(a::MultiValue{Tuple{D,D}}) where D = @notimplemented
transpose(a::MultiValue{Tuple{D,D}}) where D = @notimplemented

@generated function adjoint(a::TensorValue{D1,D2,T}) where {D1,D2,T}
  str = ""
  for i in 1:D1
    for j in 1:D2
      k = (j-1)*D1 + i
      str *= "conj(a[$i,$j]), "
    end
  end
  Meta.parse("TensorValue{D2,D1,T}($str)")
end

@generated function transpose(a::TensorValue{D1,D2,T}) where {D1,D2,T}
  str = ""
  for i in 1:D1
    for j in 1:D2
      str *= "a[$i,$j], "
    end
  end
  Meta.parse("TensorValue{D2,D1,T}($str)")
end

@inline function adjoint(a::TensorValue{D1,D2,T}) where {D1,D2,T<:Real}
  transpose(a)
end

adjoint(a::AbstractSymTensorValue) = conj(a)

@inline adjoint(a::AbstractSymTensorValue{D,T} where {D,T<:Real}) = transpose(a)

transpose(a::AbstractSymTensorValue) = a

###############################################################
# Symmetric part
###############################################################

"""
    symmetric_part(v::MultiValue{Tuple{D,D}})::AbstractSymTensorValue

Return the symmetric part of second order tensor, that is `½(v + vᵀ)`.
Return `v` if  `v isa AbstractSymTensorValue`.
"""
@generated function symmetric_part(v::MultiValue{Tuple{D,D},T}) where {D,T}
  iszero(D) && return :( zero(SymTensorValue{0,T}) )
  str = "("
  for j in 1:D
      for i in j:D
          str *= "0.5*v[$i,$j] + 0.5*v[$j,$i], "
      end
  end
  str *= ")"
  Meta.parse("SymTensorValue{D}($str)")
end

symmetric_part(v::AbstractSymTensorValue) = v

"""
    skew_symmetric_part(v::MultiValue{Tuple{D,D}})::MultiValue{Tuple{D,D}}

Return the asymmetric part of second order tensor, that is `½(v - vᵀ)`.
Return `v` if  `v isa AbstractSymTensorValue`.
"""
@generated function skew_symmetric_part(v::MultiValue{Tuple{D,D},T}) where {D,T}
  iszero(D) && return :( zero(TensorValue{0,0,T}) )
  str = "("
  for j in 1:D
      for i in 1:D
          str *= "0.5*v[$i,$j] - 0.5*v[$j,$i], "
      end
  end
  str *= ")"
  Meta.parse("TensorValue{D,D}($str)")
end

###############################################################
# diag
###############################################################

function LinearAlgebra.diag(a::MultiValue{Tuple{D,D},T}) where {D,T}
  VectorValue{D,T}((a[i,i] for i in 1:D)...)
end

###############################################################
# Broadcast
###############################################################
# TODO more cases need to be added

function Base.broadcasted(f,a::VectorValue,b::VectorValue)
  VectorValue(map(f,a.data,b.data))
end

function Base.broadcasted(f,a::TensorValue,b::TensorValue)
  TensorValue(map(f,a.data,b.data))
end

function Base.broadcasted(f,a::AbstractSymTensorValue,b::AbstractSymTensorValue)
  SymTensorValue(map(f,a.data,b.data))
end
