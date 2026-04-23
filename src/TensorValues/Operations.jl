###############################################################
# Comparison
###############################################################

(==)(a::Number, b::MultiValue) = false
(==)(a::MultiValue, b::MultiValue) = false
(==)(a::MultiValue{S}, b::MultiValue{S}) where S = a.data == b.data
(≈)(a::MultiValue,b::MultiValue;kwargs...) = ≈(get_array(a),get_array(b);kwargs...)

function (≈)(
  a::AbstractArray{<:MultiValue}, b::AbstractArray{<:MultiValue}; kwargs...)
  if size(a) != size(b); return false; end
  for (ai,bi) in zip(a,b)
    if !≈(ai,bi;kwargs...); return false; end
  end
  true
end

function isless(a::MultiValue{Tuple{L}}, b::MultiValue{Tuple{L}}) where {L}
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

isless(a::Number, b::MultiValue) = all(isless.(a, b.data))
function <=(a::Number, b::MultiValue)
  all(a .<= b.data) # default doesen't work because a==b is always false
end
isless(a::MultiValue, b::MultiValue) = @unreachable "Comparison is not defined between tensor of order greater than 1"
function <=(a::MultiValue, b::MultiValue)
  isless(a,b) || isequal(a,b)
end

###############################################################
# promotion and conversions
###############################################################

# No default promotion /conversion for tensors of different type names
promote_rule(::Type{<:MultiValue}, ::Type{<:MultiValue}) = Union{}

# But promotion and conversion between the different types of square tensors.
promote_rule(::Type{<:TensorValue{D,D,Ta}},     ::Type{<:SymTensorValue{D,Tb}})          where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}
promote_rule(::Type{<:TensorValue{D,D,Ta}},     ::Type{<:SkewSymTensorValue{D,Tb}})      where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}
promote_rule(::Type{<:TensorValue{D,D,Ta}},     ::Type{<:SymTracelessTensorValue{D,Tb}}) where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}
promote_rule(::Type{<:SymTensorValue{D,Ta}},    ::Type{<:SkewSymTensorValue{D,Tb}})      where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}
promote_rule(::Type{<:SymTensorValue{D,Ta}},    ::Type{<:SymTracelessTensorValue{D,Tb}}) where {D,Ta,Tb} = SymTensorValue{D,promote_type(Ta,Tb)}
promote_rule(::Type{<:SkewSymTensorValue{D,Ta}},::Type{<:SymTracelessTensorValue{D,Tb}}) where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}

convert(::Type{<:TensorValue{D,D,Ta}},  a::SymTensorValue{D,Tb})          where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}(get_array(a))
convert(::Type{<:TensorValue{D,D,Ta}},  a::SkewSymTensorValue{D,Tb})      where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}(get_array(a))
convert(::Type{<:TensorValue{D,D,Ta}},  a::SymTracelessTensorValue{D,Tb}) where {D,Ta,Tb} = TensorValue{D,D,promote_type(Ta,Tb)}(get_array(a))
convert(::Type{<:SymTensorValue{D,Ta}}, a::SymTracelessTensorValue{D,Tb}) where {D,Ta,Tb} = SymTensorValue{D,promote_type(Ta,Tb)}(a.data)

"""
    const _Scalar = Union{Real,Complex}

Abstract type for the scalar types that `MultiValue` support component-wise
operations with.
"""
const _Scalar = Union{Real,Complex}

# TODO deprecate next two methods ? A few stuff depend on this behavior but very ugly.
function convert(V::Type{<:MultiValue}, a::T) where T<:_Scalar
  isone(length(V)) && isone(num_indep_components(V)) || error("Cannot convert value of type $V to type $T")
  V(a)
end
function convert(T::Type{<:_Scalar}, a::V) where V<:MultiValue
  isone(length(a)) || error("Cannot convert value of type $V to type $T")
  T(a[1])
end

###############################################################
# Addition / subtraction
###############################################################

Base.iszero(a::MultiValue) = all(iszero.(a.data))

for op in (:+,:-)
  @eval begin

    function ($op)(a::T) where T<:MultiValue
      Li = num_indep_components(T)
      r = map($op, Tuple(a)[1:Li])
      T(r)
    end

    function ($op)(a::V, b::V) where V<:MultiValue
      Li = num_indep_components(V)
      r = map(($op), Tuple(a)[1:Li], Tuple(b)[1:Li])
      V(r)
    end
  end
end

###############################################################
# Matrix Division
###############################################################

function (\)(a::MultiValue{Tuple{D,D}} where D, b::MultiValue)
  r = get_array(a) \ get_array(b)
  T = change_eltype(b, eltype(r))
  T(r)
end

###############################################################
# Operations with other numbers
###############################################################

@generated function _bc(f, a::NTuple{N}, b::Number) where N
  s = "("
  for i in 1:N
    s *= "f(a[$i],b), "
  end
  s *= ")"
  Meta.parse(s)
end

@generated function _bc(f, b::Number, a::NTuple{N}) where N
  s = "("
  for i in 1:N
    s *= "f(b,a[$i]), "
  end
  s *= ")"
  Meta.parse(s)
end

for op in (:+,:-,:*)
  @eval begin
    function ($op)(a::MultiValue, b::_Scalar)
      Li = num_indep_components(a)
      r = _bc($op, Tuple(a)[1:Li], b)
      T = _eltype($op, r, a, b)
      M = change_eltype(a, T)
      M(r)
    end

    function ($op)(a::_Scalar, b::MultiValue)
      Li = num_indep_components(b)
      r = _bc($op, a, Tuple(b)[1:Li])
      T = _eltype($op, r, a, b)
      M = change_eltype(b, T)
      M(r)
    end
  end
end

const _AbstractTracelessTensor{D} = Union{SymTracelessTensorValue{D},SkewSymTensorValue{D}}
_err = "This operation is undefined for traceless tensors"
(+)(::_AbstractTracelessTensor, ::_Scalar) = error(_err)
(+)(::_Scalar, ::_AbstractTracelessTensor) = error(_err)
(-)(::_AbstractTracelessTensor, ::_Scalar) = error(_err)
(-)(::_Scalar, ::_AbstractTracelessTensor) = error(_err)

function (/)(a::MultiValue, b::_Scalar)
  Li = num_indep_components(a)
  r = _bc(/, Tuple(a)[1:Li], b)
  T = _eltype(/, r, a, b)
  P = change_eltype(a, T)
  P(r)
end

@inline function _eltype(op, r, a)
  eltype(r)
end

@inline function _eltype(op, r, a, b)
  eltype(r)
end

@inline function _eltype(op, r, a...)
  eltype(r)
end

@inline function _eltype(op, r::Tuple{}, a)
  typeof(op(zero(eltype(a))))
end

@inline function _eltype(op, r::Tuple{}, a, b)
  typeof(op(zero(eltype(a)), zero(eltype(b))))
end

@inline function _eltype(op, r::Tuple{}, a...)
  typeof(reduce(op, zero.(eltype.(a))))
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

"""
    dot(a::MultiValue{Tuple{...,D}}, b::MultiValue{Tuple{D,...}})
    a ⋅¹ b
    a ⋅ b

Single contraction of two tensors `a` and `b`, of the last index of `a` with
the first index of `b`. The corresponding dimensions `D` must match. No symmetry
is preserved.
On two vectors, this is the same as the inner product.
"""
dot(a::MultiValue, b::MultiValue) = contracted_product(Val(1), a, b)

const ⋅¹ = dot

###############################################################
# Inner product (full contraction)
###############################################################

inner(a::_Scalar, b::_Scalar) = a * b

"""
    inner(a::MultiValue{S}, b::MultiValue{S}) -> scalar
    a ⊙ b

Inner product of two tensors, that is the full contraction along each indices. The size `S` of `a` and `b` must match.
"""
function inner(a::MultiValue, b::MultiValue)
  @notimplemented "Sizes of tensors must match."
end

inner(a::MultiValue{S,Ta,N}, b::MultiValue{S,Tb,N}) where {S,Ta,Tb,N} = contracted_product(Val(N),a,b)

@generated function inner(a::AbstractSymTensorValue{D,Ta}, b::AbstractSymTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :(zero($(Base.promote_op(*,Ta,Tb))))
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

@generated function inner(a::SymFourthOrderTensorValue{D,Ta}, b::SymFourthOrderTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :(zero($(Base.promote_op(*,Ta,Tb))))

  S = Tuple{D,D,D,D}
  VInt = change_eltype(a, Int)
  # each independent component appear either 1,2 of 4 times in inputs
  strs = Dict(1 => "(", 2 => "2*(", 4 => "4*(")

  # for each independent component, add its product in the sum of the
  # corresponding multiplicative factor
  for (i, fi) in enumerate(component_basis(VInt))
    factor = @invoke inner(fi::MultiValue{S,Int,4}, fi::MultiValue{S,Int,4}) # use the generic but slower method
    strs[factor] *= "+ indep_comp_getindex(a,$i)*indep_comp_getindex(b,$i)"
  end

  str = string(strs[1][:], ") + ", strs[2][:], ") + ", strs[4][:], ")")
  Meta.parse(str)
end

function inner(a::SkewSymTensorValue{D,Ta}, b::SkewSymTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return zero(Base.promote_op(*,Ta,Tb))
  2 * inner(VectorValue(a.data), VectorValue(b.data))
end

function inner(a::SkewSymTensorValue{D,Ta}, b::AbstractSymTensorValue{D,Tb}) where {D,Ta,Tb}
  zero(Base.promote_op(*,Ta,Tb))
end
function inner(a::AbstractSymTensorValue{D,Tb}, b::SkewSymTensorValue{D,Ta}) where {D,Ta,Tb}
  zero(Base.promote_op(*,Ta,Tb))
end

# TODO These two methods make no sense and shold be removed
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
double_contraction(a::MultiValue, b::MultiValue) = contracted_product(Val(2), a, b)

# c_ijkl = a_ijmn*b_mnkl (3D)
@generated function double_contraction(a::SymFourthOrderTensorValue{3}, b::SymFourthOrderTensorValue{3})

  Sym4TensorIndexing = [1111, 1121, 1131, 1122, 1132, 1133, 2111, 2121, 2131, 2122, 2132, 2133,
    3111, 3121, 3131, 3122, 3132, 3133, 2211, 2221, 2231, 2222, 2232, 2233,
    2311, 2321, 2331, 2322, 2332, 2333, 3311, 3321, 3331, 3322, 3332, 3333]
  ss = String[]
  for off_index in Sym4TensorIndexing
    i = parse(Int, string(off_index)[1])
    j = parse(Int, string(off_index)[2])
    k = parse(Int, string(off_index)[3])
    l = parse(Int, string(off_index)[4])
    s = join(["a[$j,$i,$m,$n]*b[$m,$n,$l,$k]+" for m in 1:3 for n in 1:3])
    push!(ss, s[1:(end-1)] * ", ")
  end
  str = join(ss)
  Meta.parse("SymFourthOrderTensorValue{3}($str)")
end

function _comp_prod_double_symfourth4(::Val{D},a,b,i,j,k,l) where D
  s = zero(promote_type(eltype(a), eltype(b)))
  @inbounds for m in 1:D
    for n in 1:D
      s += a[i,j,m,n]*b[m,n,k,l]
    end
  end
  s
end

# c_ijkl = a_ijmn*b_mnkl (general case)
@generated function double_contraction(a::SymFourthOrderTensorValue{D,Ta}, b::SymFourthOrderTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :(SymFourthOrderTensorValue{0,$(Base.promote_op(*,Ta,Tb))}())

  str = ""
  for j in 1:D
    for i in j:D
      for l in 1:D
        for k in l:D
          if D < 4
            s = ""
            for m in 1:D
              for n in 1:D
                s *= " a[$i,$j,$m,$n]*b[$m,$n,$k,$l] +"
              end
            end
            str *= s[1:(end-1)] * ", "
          else
            # for D=4, this compiles in <0.05 s, while the other takes 200 s (Julia 1.12).
            # But runtime is 5 times slower 1.6 μs, vs 300 ns
            #
            # I tried passing i,j,k,l by Val() instead of value, compilation
            # takes 1s, the runtime isn't improved
            #
            # This means that the acceleration we get by removing the function
            # likely comes from llvm optimization that recycle partial results
            # of the products for different i,j,k,l
            str *= "_comp_prod_double_symfourth4(Val(D),a,b,$i,$j,$k,$l), "
          end
        end
      end
    end
  end
  Meta.parse("SymFourthOrderTensorValue{D}($str)")
end

# c_ijk = a_ilm*b_lmjk
@generated function double_contraction(a::ThirdOrderTensorValue{D1,D,D,Ta}, b::SymFourthOrderTensorValue{D,Tb}) where {D1,D,Ta,Tb}
  iszero(length(a)) && return :(zero(ThirdOrderTensorValue{D1,D,D,$(Base.promote_op(*,Ta,Tb))}))
  ss = String[]
  for k in 1:D
    for j in 1:D
      for i in 1:D1
        s = join(["a[$i,$l,$m]*b[$l,$m,$j,$k]+" for l in 1:D for m in 1:D])
        push!(ss, s[1:(end-1)] * ", ")
      end
    end
  end
  str = join(ss)
  Meta.parse("ThirdOrderTensorValue{$D1,$D,$D}($str)")
end

# c_ij = a_ijkl*b_kl
@generated function double_contraction(a::SymFourthOrderTensorValue{D,Ta}, b::AbstractSymTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :(zero(SymTensorValue{D,$(Base.promote_op(*,Ta,Tb))}))
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

# c_ij = a_kl*b_klij
@generated function double_contraction(a::AbstractSymTensorValue{D,Ta}, b::SymFourthOrderTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :(zero(SymTensorValue{D,$(Base.promote_op(*,Ta,Tb))}))
  str = ""
  for i in 1:D
    for j in i:D
      for k in 1:D
        str *= "+ a[$k,$k]*b[$k,$k,$i,$j]"
      end
      str *= " + 2*("
      for k in 1:D
        for l in k+1:D
          str *= "+ a[$k,$l]*b[$k,$l,$i,$j]"
        end
      end
      str *= "), "
    end
  end
  Meta.parse("SymTensorValue{D}($str)")
end

# c_ij = a_ijkl*b_kl
function double_contraction(a::SymFourthOrderTensorValue{D}, b::MultiValue{Tuple{D,D}}) where {D}
  double_contraction(a, symmetric_part(b))
end

# c_ij = a_kl*b_klij
function double_contraction(a::MultiValue{Tuple{D,D}}, b::SymFourthOrderTensorValue{D}) where {D}
  double_contraction(symmetric_part(a), b)
end

const ⋅² = double_contraction

###############################################################
# Congruent product
###############################################################

"""
    congruent_prod(a, b)

Given a square second order tensors `a` and `b`, return `b`ᵀ⋅`a`⋅`b`.
The type of the resulting value is (skew) symmetric stable w.r.t. `typeof(a)`.
"""
function congruent_prod(a::MultiValue{Tuple{D,D},Ta}, b::MultiValue{Tuple{D,D1},Tb}) where {D,D1,Ta,Tb}
  T = Base.promote_op(*,Ta,Tb)
  V = _congruent_ret_type(a, D1)
  (iszero(D) || iszero(D1)) && return zero(V{T})
  V(get_array(transpose(b) ⋅ a ⋅ b))
end
_congruent_ret_type(a, D1) = TensorValue{D1,D1}
_congruent_ret_type(a::AbstractSymTensorValue, D1) = SymTensorValue{D1}
_congruent_ret_type(a::SkewSymTensorValue, D1) = SkewSymTensorValue{D1}

function congruent_prod(a::Number, b::Number)
  msg = """ operation only defined for 2nd order tensors `a` and `b` with
      `size(b,1) == size(a, 1) == size(a, 2)`, got `size(a)=$(size(a))` and
      `size(b)=$(size(b))`.
      """
  @unreachable msg
end

###############################################################
# Reductions
###############################################################

for op in (:sum, :maximum, :minimum)
  @eval begin
    $op(a::MultiValue) = $op(a.data)
  end
end

###############################################################
# Outer product (aka tensor product, also dyadic product)
###############################################################

outer(a::Number, b::Number) = a * b

outer(a::MultiValue, b::Number) = a * b
outer(a::Number, b::MultiValue) = a * b

"""
    outer(a,b)
    a ⊗ b

Outer product (or tensor-product) of two `Number`s and/or `MultiValue`s, that is
`(a⊗b)[i₁,...,iₙ,j₁,...,jₙ] = a[i₁,...,iₙ]*b[j₁,...,jₙ]`. This falls back to standard
multiplication if `a` or `b` is a scalar.
"""
outer(a::MultiValue, b::MultiValue) = contracted_product(Val(0), a, b)

# c_ijkl = a_ij*b_kl
@generated function outer(a::AbstractSymTensorValue{D,Ta}, b::AbstractSymTensorValue{D,Tb}) where {D,Ta,Tb}
  iszero(D) && return :(zero(SymFourthOrderTensorValue{D,$(Base.promote_op(*,Ta,Tb))}))
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
# Generic tensor product + contraction
###############################################################

"""
    contracted_product(Val(n), a::MultiValue, b::MultiValue)

Compute the tensor `c` that results from the tensor-product of `a` and `b`
contracted over the `n` last indices of `a` and `n` first indices of `b`:

`c[i₁,...,iₘ₋ₙ,jₙ₊₁,...,jₚ] = Σ_{l₁,...,lₙ}  a[i₁,...,iₘ₋ₙ,l₁,...,lₙ]*b[l₁,...,lₙ,jₙ₊₁,...,jₚ]`

where `m` is the order `a` and `p` that of `b`. The sum runs for each cartesian
indices `(l₁,...,lₙ)` of the last `n` dimension of `a`, which must match the
first `n` dimensions of `b`.

This operation generalizes:
- [`outer` / `⊗`](@ref outer) for `n=0`,
- [`dot` / `⋅`](@ref dot) for `n=1`,
- [`double_contraction` / `⋅²`](@ref double_contraction) for `n=2`,
- [`inner` / `⊙`](@ref inner) for `n = length(size(a)) = length(size(b))`.

This function is not optimized for tensor types with dependent components (use
the specific functions above if possible), but is used as default generic implementation.
"""
@generated function contracted_product(
  ::Val{N}, a::MultiValue{Sa,Ta,Na}, b::MultiValue{Sb,Tb,Nb}) where {N,Sa,Ta,Na,Sb,Tb,Nb}

  @assert ( Sa <: Tuple && Sb <: Tuple ) "Ill-defined MultiValue value"

  # check that no tensor is of order 0.
  if min(Na,Nb) == 0
    msg =  "Generic $N-contraction only implemented if both tensors are of order ≥ 1,
      got tensor sizes $Sa and $Sb"
    return :( error($msg) )
  end

  last_N_Sa = last(Sa.parameters, N)
  first_N_Sb = first(Sb.parameters, N)
  if ( length(last_N_Sa)≠N || length(first_N_Sb)≠N || last_N_Sa != first_N_Sb )
    msg = """$N-contraction is not possible for tensors of size $Sa and $Sb, it is required that
      - both tensor orders (size lengths) must be greater than or equal to $N
      - the last $N dimension of the first argument match the first $N of the second argument
      """
    return :( throw(DimensionMismatch($(msg) )) )
  end

  S_contract = tuple(last_N_Sa...)
  Sa_keep, Sb_keep = tuple(Sa.parameters[1:Na-N]...), tuple(Sb.parameters[N+1:Nb]...)

  Sr = tuple(Sa_keep..., Sb_keep...)
  Nr = length(Sr)
  Vstr = if Nr == 0
    ""
  elseif Nr == 1
    "VectorValue{$(Sr[1])}"
  elseif Nr == 2
    "TensorValue{$(Sr[1]),$(Sr[2])}"
  else
    "HighOrderTensorValue{$(Tuple{Sr...})}"
  end

  if (iszero(length(a)) || iszero(length(b)))
    if iszero(Nr)
      return Meta.parse("zero(Base.promote_op(*,Ta,Tb))")
    else
      return Meta.parse("zero("*Vstr*"{Base.promote_op(*,Ta,Tb)})")
    end
  end

  # TODO, if the length of the resulting tensor exceeds some threshold, switch
  # to runtime looping over the indices, or even using BLAS. Also, we might need
  # to change HighOrderTensorValue to store into Memory or simply some fixed sized array.
  # Context: compiling double_contaction of 4th order 4D tensors takes 5-10 min, runtime 15μs.
  ss = String[Vstr, "("]
  for cib in CartesianIndices(Sb_keep) # Julia is column major, last index enumerates first
    for cia in CartesianIndices(Sa_keep)
      push!(ss, join("+a[$cia, $ciC]*b[$ciC, $cib]" for ciC in CartesianIndices(S_contract)))
      push!(ss, ", ")
    end
  end
  pop!(ss) #rm last comma in case of scalar output
  push!(ss, ")")
  Meta.parse(join(ss))
end

###############################################################
# General tensor contraction (arbitrary index pairs)
###############################################################

"""
    tensor_contraction(a::MultiValue, b::MultiValue, ia::NTuple{N,Int}, ib::NTuple{N,Int})

Contract `N` index pairs between tensors `a` and `b`: the `ia[k]`-th index of `a`
is summed against the `ib[k]`-th index of `b`, for `k = 1, …, N`.
The contracted dimensions must match pairwise.

The output tensor has the remaining (non-contracted) indices of `a` in their
original order, followed by those of `b` in their original order.

Generalises [`contracted_product`](@ref): `contracted_product(Val(N), a, b)` is
equivalent to `tensor_contraction(a, b, (Na-N+1, …, Na), (1, …, N))`.
"""
function tensor_contraction(
  a::MultiValue{Sa,Ta,Na},
  b::MultiValue{Sb,Tb,Nb},
  ia::NTuple{Nc,Int},
  ib::NTuple{Nc,Int}) where {Sa,Ta,Na,Sb,Tb,Nb,Nc}
  _tensor_contraction(a, b, Val(ia), Val(ib))
end

tensor_contraction(a::MultiValue, b::MultiValue, ia::Int, ib::Int) =
  tensor_contraction(a, b, (ia,), (ib,))

@generated function _tensor_contraction(
  a::MultiValue{Sa,Ta,Na},
  b::MultiValue{Sb,Tb,Nb},
  ::Val{Ia},
  ::Val{Ib}) where {Sa,Ta,Na,Sb,Tb,Nb,Ia,Ib}

  @assert (Sa <: Tuple && Sb <: Tuple) "Ill-defined MultiValue"
  Nc = length(Ia)
  @assert length(Ib) == Nc

  if min(Na, Nb) == 0
    msg = "tensor_contraction requires tensors of order ≥ 1, got orders $Na and $Nb"
    return :(error($msg))
  end

  for i in Ia
    (1 <= i <= Na) || begin
      msg = "ia index $i out of range [1, $Na]"
      return :(throw(ArgumentError($msg)))
    end
  end
  for i in Ib
    (1 <= i <= Nb) || begin
      msg = "ib index $i out of range [1, $Nb]"
      return :(throw(ArgumentError($msg)))
    end
  end
  length(unique(Ia)) == Nc || return :(throw(ArgumentError("ia contains duplicate indices")))
  length(unique(Ib)) == Nc || return :(throw(ArgumentError("ib contains duplicate indices")))

  for k in 1:Nc
    da, db = Sa.parameters[Ia[k]], Sb.parameters[Ib[k]]
    if da != db
      msg = "Dimension mismatch at contracted pair $k: a[$(Ia[k])]=$da, b[$(Ib[k])]=$db"
      return :(throw(DimensionMismatch($msg)))
    end
  end

  Ia_set = Set(Ia)
  Ib_set = Set(Ib)
  ka = tuple(filter(i -> !(i in Ia_set), 1:Na)...)
  kb = tuple(filter(i -> !(i in Ib_set), 1:Nb)...)

  S_a_kept   = ntuple(j -> Sa.parameters[ka[j]], length(ka))
  S_b_kept   = ntuple(j -> Sb.parameters[kb[j]], length(kb))
  S_contract = ntuple(k -> Sa.parameters[Ia[k]], Nc)

  Sr = (S_a_kept..., S_b_kept...)
  Nr = length(Sr)

  Vstr = if Nr == 0
    "Base.promote_op(*,Ta,Tb)"
  elseif Nr == 1
    "VectorValue{$(Sr[1]),Base.promote_op(*,Ta,Tb)}"
  elseif Nr == 2
    "TensorValue{$(Sr[1]),$(Sr[2]),Base.promote_op(*,Ta,Tb)}"
  else
    "HighOrderTensorValue{$(Tuple{Sr...}),Base.promote_op(*,Ta,Tb)}"
  end

  if iszero(length(a)) || iszero(length(b))
    return Meta.parse("zero($Vstr)")
  end

  s_a_ranges = ntuple(j -> 1:S_a_kept[j], length(ka))
  s_b_ranges = ntuple(j -> 1:S_b_kept[j], length(kb))
  s_c_ranges = ntuple(k -> 1:S_contract[k], Nc)

  ss = String[Vstr * "("]
  for cib in Iterators.product(s_b_ranges...)
    for cia in Iterators.product(s_a_ranges...)
      terms = String[]
      for ciC in Iterators.product(s_c_ranges...)
        a_idx = Vector{Int}(undef, Na)
        for (j, pos) in enumerate(ka)
          a_idx[pos] = cia[j]
        end
        for (k, pos) in enumerate(Ia)
          a_idx[pos] = ciC[k]
        end
        b_idx = Vector{Int}(undef, Nb)
        for (j, pos) in enumerate(kb)
          b_idx[pos] = cib[j]
        end
        for (k, pos) in enumerate(Ib)
          b_idx[pos] = ciC[k]
        end
        push!(terms, "+a[$(join(a_idx,","))]*b[$(join(b_idx,","))]")
      end
      push!(ss, join(terms) * ", ")
    end
  end
  push!(ss, ")")
  Meta.parse(join(ss))
end

"""
    tensor_contraction(a::MultiValue, i::NTuple{N,Int}, j::NTuple{N,Int})

Self-contraction of `a`: for each pair `k`, index `i[k]` and index `j[k]` of `a`
are set equal and summed over. `i` and `j` must be disjoint and the paired
dimensions must match.

The output has the remaining indices of `a` (those not in `i` or `j`) in their
original order. `tr(a)` is the special case `tensor_contraction(a, (1,), (2,))`.
"""
function tensor_contraction(
  a::MultiValue{Sa,Ta,Na},
  i::NTuple{Nc,Int},
  j::NTuple{Nc,Int}) where {Sa,Ta,Na,Nc}
  _tensor_contraction(a, Val(i), Val(j))
end

tensor_contraction(a::MultiValue, i::Int, j::Int) =
  tensor_contraction(a, (i,), (j,))

@generated function _tensor_contraction(
  a::MultiValue{Sa,Ta,Na},
  ::Val{I},
  ::Val{J}) where {Sa,Ta,Na,I,J}

  @assert Sa <: Tuple "Ill-defined MultiValue"
  Nc = length(I)
  @assert length(J) == Nc

  Na == 0 && return :(error("tensor_contraction requires a tensor of order ≥ 1"))

  for idx in I
    (1 <= idx <= Na) || begin
      msg = "i index $idx out of range [1, $Na]"
      return :(throw(ArgumentError($msg)))
    end
  end
  for idx in J
    (1 <= idx <= Na) || begin
      msg = "j index $idx out of range [1, $Na]"
      return :(throw(ArgumentError($msg)))
    end
  end
  length(unique([I..., J...])) == 2Nc ||
    return :(throw(ArgumentError("i and j must be disjoint with no internal duplicates")))

  for k in 1:Nc
    di, dj = Sa.parameters[I[k]], Sa.parameters[J[k]]
    di == dj || begin
      msg = "Dimension mismatch at pair $k: a[$(I[k])]=$di ≠ a[$(J[k])]=$dj"
      return :(throw(DimensionMismatch($msg)))
    end
  end

  contracted = Set([I..., J...])
  ka         = tuple(filter(idx -> !(idx in contracted), 1:Na)...)
  S_keep     = ntuple(m -> Sa.parameters[ka[m]], length(ka))
  S_contract = ntuple(k -> Sa.parameters[I[k]], Nc)
  Nr         = length(S_keep)

  Vstr = if Nr == 0
    "Ta"
  elseif Nr == 1
    "VectorValue{$(S_keep[1]),Ta}"
  elseif Nr == 2
    "TensorValue{$(S_keep[1]),$(S_keep[2]),Ta}"
  else
    "HighOrderTensorValue{$(Tuple{S_keep...}),Ta}"
  end

  iszero(length(a)) && return Meta.parse("zero($Vstr)")

  s_keep_ranges     = ntuple(m -> 1:S_keep[m],     length(ka))
  s_contract_ranges = ntuple(k -> 1:S_contract[k], Nc)

  ss = String[Vstr * "("]
  for cia in Iterators.product(s_keep_ranges...)
    terms = String[]
    for ciC in Iterators.product(s_contract_ranges...)
      a_idx = Vector{Int}(undef, Na)
      for (m, pos) in enumerate(ka)
        a_idx[pos] = cia[m]
      end
      for (k, pos) in enumerate(I)
        a_idx[pos] = ciC[k]
      end
      for (k, pos) in enumerate(J)
        a_idx[pos] = ciC[k]
      end
      push!(terms, "+a[$(join(a_idx, ","))]")
    end
    push!(ss, join(terms) * ", ")
  end
  push!(ss, ")")
  Meta.parse(join(ss))
end

###############################################################
# Cross Product
###############################################################

function cross(a::MultiValue{Tuple{3}}, b::MultiValue{Tuple{3}})
  VectorValue{3}(a[2]b[3] - a[3]b[2], a[3]b[1] - a[1]b[3], a[1]b[2] - a[2]b[1])
end

function cross(a::MultiValue{Tuple{2}}, b::MultiValue{Tuple{2}})
  a[1]b[2] - a[2]b[1]
end

"""
    cross(a::VectorValue{3}, b::VectorValue{3}) -> VectorValue{3}
    cross(a::VectorValue{2}, b::VectorValue{2}) -> scalar
    a × b

Cross product of 2D and 3D vector.
"""
cross(a::MultiValue, b::MultiValue) = error("Cross product only defined for R2 and R3 vectors of same dimension")

###############################################################
# Linear Algebra
###############################################################

"""
    det(a::MultiValue{Tuple{D,D},T})

Determinent of square second order tensors.
"""
det(a::MultiValue{Tuple{D,D}}) where {D} = det(get_array(a))
det(a::MultiValue) = @unreachable "det undefined for this tensor shape: $(size(a))"

det(a::MultiValue{Tuple{1,1}}) = a[1]

function det(a::MultiValue{Tuple{2,2}})
  a_11 = a[1, 1]; a_12 = a[1, 2]
  a_21 = a[2, 1]; a_22 = a[2, 2]
  a_11*a_22 - a_12*a_21
end

function det(a::MultiValue{Tuple{3,3}})
  a_11 = a[1,1]; a_12 = a[1,2]; a_13 = a[1,3]
  a_21 = a[2,1]; a_22 = a[2,2]; a_23 = a[2,3]
  a_31 = a[3,1]; a_32 = a[3,2]; a_33 = a[3,3]
  a_11*a_22*a_33 + a_12*a_23*a_31 + a_13*a_21*a_32 -
    (a_11*a_23*a_32 + a_12*a_21*a_33 + a_13*a_22*a_31)
end

det(::SkewSymTensorValue{3,T}) where T = zero(T)

"""
    inv(a::MultiValue{Tuple{D,D}})

Inverse of a second order tensor.
"""
inv(a::MultiValue{Tuple{D,D}}) where D = TensorValue(inv(get_array(a)))

const InverseStableTensorTypes{D} = Union{SymTensorValue{D},SkewSymTensorValue{D}}

function inv(a::InverseStableTensorTypes{D}) where D
  ai = inv(get_array(a))
  T = change_eltype(a, eltype(ai))
  T(ai)
end

function inv(a::MultiValue{Tuple{1,1}})
  r = 1 / a[1]
  T = change_eltype(a, typeof(r))
  T(r)
end

function inv(a::MultiValue{Tuple{2,2}})
  c = 1 / det(a)
  data = (a[2, 2] * c, -a[2, 1] * c, -a[1, 2] * c, a[1, 1] * c)
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

function inv(a::SymTensorValue{2})
  c = 1/det(a)
  T = change_eltype(a,typeof(c))
  T(a[2,2]*c, -a[2,1]*c, a[1,1]*c)
end

inv(::SymTracelessTensorValue{1,T}) where T = TensorValue{1,1}(inv(zero(T)))
function inv(a::SymTracelessTensorValue{2})
  c = -1/det(a)
  T = change_eltype(a,typeof(c))
  T(a[1,1]*c, a[2,1]*c)
end

inv(::SkewSymTensorValue{1,T}) where T = TensorValue{1,1}(inv(zero(T)))
inv(a::SkewSymTensorValue{2}) = (typeof(a))(-inv(a.data[1]))
inv(a::SkewSymTensorValue{3,T,L}) where {T,L} = SkewSymTensorValue{3,T}(tfill(inv(zero(T)), Val(L)))

"""
    eigen(a::MultiValue{Tuple{D,D}})

Eigenvalue decomposition of a square second order tensor.
"""
eigen(a::MultiValue{Tuple{D,D}}) where D = eigen(get_array(a))
eigen(a::MultiValue) = @unreachable "eigen undefined for this tensor shape: $(size(a))"

"""
    normalize(a::MultiValue)

Normalization of a tensor value.
"""
normalize(a::MultiValue) = a / norm(a)

sqrt(t::MultiValue{Tuple{D,D}}) where D = TensorValue{D,D}(sqrt(SArray(t)))

###############################################################
# Measure
###############################################################

"""
    meas(a::MultiValue{Tuple{D}})
    meas(a::MultiValue{Tuple{1,D2}})

Euclidean norm of a vector.
"""
meas(a::MultiValue{Tuple{D}}) where D = sqrt(inner(a, a))

"""
    meas(J::MultiValue{Tuple{D1,D2}})

Returns the absolute `D1`-dimensional volume of the parallelepiped
formed by the rows of `J`, that is `sqrt(det(J⋅Jᵀ))`, or `abs(det(J))` if `D1`=`D2`.
This is used to compute the contribution of the Jacobian matrix `J` of a changes of variables in integrals.
"""
meas(a::MultiValue{Tuple{D,D}}) where D = abs(det(a))

function meas(v::MultiValue{Tuple{1,D}}) where D
  t = VectorValue(v.data)
  meas(t)
end

function meas(v::MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1, n2, n3)
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
@inline norm(u::MultiValue{Tuple{D},<:Real}) where D = sqrt(inner(u, u))
@inline norm(u::MultiValue{Tuple{D}}) where D = sqrt(real(inner(u, conj(u))))
@inline norm(u::MultiValue{Tuple{D1,D2},<:Real}) where {D1,D2} = sqrt(inner(u, u))
@inline norm(u::MultiValue{Tuple{D1,D2}}) where {D1,D2} = sqrt(real(inner(u, conj(u))))
@inline norm(u::MultiValue{Tuple{0},T}) where T<:Real= sqrt(zero(T))
@inline norm(u::MultiValue{Tuple{0},T}) where T = sqrt(real(zero(T)))

###############################################################
# conj, real, imag
###############################################################

for op in (:conj,:real,:imag)
  @eval begin
    function ($op)(a::T) where T<:MultiValue
      Li = num_indep_components(a)
      r = map($op, Tuple(a)[1:Li])
      T2 = _eltype($op,r,a)
      M = change_eltype(a,T2)
      M(r)
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
  str = join([" v[$i,$i] +" for i in 1:D])
  Meta.parse(str[1:(end-1)])
end
tr(::SymTracelessTensorValue{D,T}) where {D,T} = zero(T)
tr(::SkewSymTensorValue{D,T}) where {D,T} = zero(T)
tr(::MultiValue{Tuple{A,B}}) where {A,B} = throw(ArgumentError("Second order tensor is not square"))

"""
    tr(v::MultiValue{Tuple{D,D,D2}}) -> ::VectorValue{D2}

Return a vector of length `D2` of traces computed on the first two indices: `resⱼ = Σᵢ vᵢᵢⱼ`.
"""
@generated function tr(v::HighOrderTensorValue{S,T,N}) where {S,T,N}
  @assert S <: Tuple "Ill-defined MultiValue value"

  if N > 3
    # This may be generalizable if wanted
    msg = "Trace is not implemented for tensors of order ≥ 4, got tensors ::$(v)"
    return :( throw(ArgumentError($msg)) )
  end

  A12 = first(S.parameters, 2)
  if !allequal(A12)
    msg = "First two dimensions ($A12) are not iddentical."
    return :( throw(ArgumentError($msg)) )
  end

  B = last(S.parameters)
  iszero(length(v)) && return :(zero(VectorValue{$B,T}))
  str = ""
  for k in 1:B
    for i in 1:A12[1]
      if i != 1
        str *= " + "
      end
      str *= " v[$i,$i,$k]"
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

@generated function adjoint(a::TensorValue{D1,D2,T}) where {D1,D2,T}
  str = ""
  for i in 1:D1
    for j in 1:D2
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
adjoint(a::SkewSymTensorValue) = -conj(a)

transpose(a::AbstractSymTensorValue) = a
transpose(a::SkewSymTensorValue) = -a

###############################################################
# permutedims
###############################################################

"""
    permutedims(a::MultiValue{Sa,Ta,Na}, perm::NTuple{Na,Int})

Return a tensor whose indices are those of `a` reordered by `perm`:
`result[i₁,…,iₙ] = a[i_{σ⁻¹(1)},…,i_{σ⁻¹(n)}]`, where `σ = perm`.

The output shape is `(size(a, perm[1]), …, size(a, perm[N]))`.
For second-order tensors this is equivalent to `transpose` (with `perm = (2,1)`).
Symmetry of the input tensor is not preserved in the output type.
"""
function Base.permutedims(
  a::MultiValue{Sa,Ta,Na},
  perm::NTuple{Na,Int}) where {Sa,Ta,Na}
  _permutedims(a, Val(perm))
end

@generated function _permutedims(
  a::MultiValue{Sa,Ta,Na}, ::Val{P}) where {Sa,Ta,Na,P}

  @assert Sa <: Tuple "Ill-defined MultiValue"

  sort([P...]) == collect(1:Na) || begin
    msg = "$P is not a valid permutation of 1:$Na"
    return :(throw(ArgumentError($msg)))
  end

  Sa_params = tuple(Sa.parameters...)
  Sr = ntuple(k -> Sa_params[P[k]], Na)

  Vstr = if Na == 1
    "VectorValue{$(Sr[1]),Ta}"
  elseif Na == 2
    "TensorValue{$(Sr[1]),$(Sr[2]),Ta}"
  else
    "HighOrderTensorValue{$(Tuple{Sr...}),Ta}"
  end

  iszero(length(a)) && return Meta.parse("zero($Vstr)")

  inv_perm = zeros(Int, Na)
  for k in 1:Na
    inv_perm[P[k]] = k
  end

  ss = String[Vstr * "("]
  for ci in CartesianIndices(Sr)
    a_idx = ntuple(k -> ci[inv_perm[k]], Na)
    push!(ss, "a[$(join(a_idx, ","))], ")
  end
  push!(ss, ")")
  Meta.parse(join(ss))
end

###############################################################
# Symmetric and Skew symmetric parts
###############################################################

"""
    symmetric_part(v::MultiValue{Tuple{D,D}})::AbstractSymTensorValue

Return the symmetric part of second order tensor, that is `½(v + vᵀ)`.
Return `v` if  `v isa AbstractSymTensorValue`, and the zero symmetric tensor if
`v  isa SkewSymTensorValue`.
"""
@generated function symmetric_part(v::MultiValue{Tuple{D,D},T}) where {D,T}
  iszero(D) && return :(zero(SymTensorValue{0,T}))
  str = "("
  for j in 1:D
    for i in j:D
      str *= "(v[$i,$j] + v[$j,$i])/2, "
    end
  end
  str *= ")"
  Meta.parse("SymTensorValue{D}($str)")
end

symmetric_part(v::AbstractSymTensorValue) = v
symmetric_part(::SkewSymTensorValue{D,T}) where {D,T} = zero(SymTensorValue{D,T})

"""
    skew_symmetric_part(v::MultiValue{Tuple{D,D}})::SkewSymTensorValue{D}

Return the asymmetric part of `v`, that is `½(v - vᵀ)`.
Return the zero skew symmetric tensor if `v  isa AbstractSymTensorValue`, and
`v` itself if  `v  isa SkewSymTensorValue`.
"""
@generated function skew_symmetric_part(v::MultiValue{Tuple{D,D},T}) where {D,T}
  iszero(D) && return :(zero(SkewSymTensorValue{0,T}))
  str = "("
  for i in 1:D
    for j in i+1:D
      str *= "(v[$i,$j] - v[$j,$i])/2, "
    end
  end
  str *= ")"
  Meta.parse("SkewSymTensorValue{D}($str)")
end

skew_symmetric_part(::AbstractSymTensorValue{D,T}) where {D,T} = zero(SkewSymTensorValue{D,T})
skew_symmetric_part(v::SkewSymTensorValue) = v

###############################################################
# diag
###############################################################

function LinearAlgebra.diag(a::MultiValue{Tuple{D,D},T}) where {D,T}
  VectorValue{D,T}((a[i,i] for i in 1:D)...)
end

function LinearAlgebra.diag(a::SkewSymTensorValue{D,T}) where {D,T}
  zero(VectorValue{D,T})
end

###############################################################
# Broadcast
###############################################################

function Base.broadcasted(f, a::VectorValue, b::VectorValue)
  VectorValue(map(f, a.data, b.data))
end

function Base.broadcasted(f, a::TensorValue, b::TensorValue)
  TensorValue(map(f, a.data, b.data))
end

function Base.broadcasted(f, a::AbstractSymTensorValue, b::AbstractSymTensorValue)
  SymTensorValue(map(f, a.data, b.data))
end
