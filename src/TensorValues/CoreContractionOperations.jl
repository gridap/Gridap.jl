# Core tensor algebra: the generic contraction and congruence kernels
#
# The operations here are the most abstract, general and hard to read amongst all:
# - congruent_prod
# - contracted_product
# - tensor_contraction

###############################################################
# Generate-time expression builders
###############################################################

# `sgn * expr`, with the sign folded in rather than emitted as a factor.
_signed(sgn::Int, expr) = sgn > 0 ? expr : :(-$expr)

# Sum of `terms`. `Expr(:call, :+)` with no term is `+()`, which errors, so an
# empty sum is emitted as an explicit `zero(T)`.
_sum_expr(terms, T) = isempty(terms) ? :(zero($T)) : Expr(:call, :+, terms...)

# Leibniz expansion of the determinant of the `K`×`K` matrix whose `(a,b)` entry
# is the expression `entry(a,b)`. `K ≥ 1`.
function _det_expr(entry, K::Int)
  terms = (_signed(levicivita(p), Expr(:call, :*, (entry(a, p[a]) for a in 1:K)...))
           for p in permutations(1:K))
  Expr(:call, :+, terms...)
end

# Determinant of the `K`×`K` submatrix of the second order tensor named by `A`,
# with rows `I` and columns `J`.
#
# The expansion has `K!` terms per call site. Past `K = 4` that makes compilation
# the dominant cost, so the emitted code calls `_minor_k` and loops instead.
function _minor_expr(A::Symbol, I, J)
  K = length(I)
  K == 0 && return :(one(eltype($A)))
  K > 4  && return :(_minor_k($A, $(Tuple(I)), $(Tuple(J))))
  _det_expr((a, b) -> :($A[$(I[a]), $(J[b])]), K)
end

# Determinant of the K×K submatrix of `A` with rows `I` and columns `J`.
function _minor_k(A, I::NTuple{K,Int}, J::NTuple{K,Int}) where K
  K == 0 && return one(eltype(A))
  det(SMatrix{K,K}(ntuple(n -> A[I[(n-1)%K+1], J[(n-1)÷K+1]], Val(K*K))))
end

# Contracts an alternating tensor's components with the `K`-th compound matrix of M:
#
#   (Mᴷ src)_I = Σ_J det(M[I,J]) src_J
#
# one sum per sorted combination `I ∈ c_out`. If `transposed`,  `J` indices M's
# rows and `I` its columns
function _alt_compound_comps(M::Symbol, c_out, c_in, T, src::Symbol, transposed::Bool=false)
  [ _sum_expr([ :( $(transposed ? _minor_expr(M, J, I) : _minor_expr(M, I, J)) *
                   indep_comp_getindex($src, $l) )
                for (l, J) in enumerate(c_in) ], T)
    for I in c_out ]
end

###############################################################
# Canonical tensor types
###############################################################

#The canonical type (with tuple dimension set) of a tensor of order `Nr` and
#size `Tuple` `Sr`. At order 1 the variance of `Sr[1]` picks the type
#constructor, the three being distinct types.
function _tensor_base_type(Nr,Sr)
  if Nr == 0
    ScalarValue
  elseif Nr == 1
    D = Sr[1]
    if D isa Int
      VectorValue{D}
    elseif D isa Cov
      CoVectorValue{D}
    elseif D isa Con
      ContraVectorValue{D}
    else
      @unreachable "$D is not a tensor dimension, see `Dim`."
    end
  elseif Nr == 2
    TensorValue{Sr[1],Sr[2]}
  else
    NTensorValue{Tuple{Sr...}}
  end
end

# expr of _tensor_base_type, except for Nr=0, empty expr is returned instead of ScalarValue
_tensor_base_type_expr(Nr,Sr) =
  Nr == 0 ? :() : :( $(_tensor_base_type(Nr,Sr)) )

###############################################################
# Generic tensor product + contraction
###############################################################

function _zero_tensor_expr(V,Nr,Ta,Tb)
    if iszero(Nr)
      :( zero(Base.promote_op(*,Ta,Tb)) )
    else
      :( zero($V{Base.promote_op(*,Ta,Tb)}) )
    end
end

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

  # An order-0 tensor can only be contracted if has exactly one component
  if (iszero(Na) && length(a) != 1) || (iszero(Nb) && length(b) != 1)
    msg = "Generic $N-contraction requires an order-0 tensor to hold exactly one
      component, got tensor sizes $Sa and $Sb"
    return :( error($msg) )
  end

  sz_a, sz_b = size(a), size(b)
  last_N_Sa = last(sz_a, N)
  first_N_Sb = first(sz_b, N)
  if ( length(last_N_Sa)≠N || length(first_N_Sb)≠N)
    msg = "$N-contraction is not possible for tensors of size $Sa and $Sb, one tensor does not have $N indices"
    return :( throw(DimensionMismatch($(msg) )) )
  end
  for (i,j) in zip(last_N_Sa, first_N_Sb)
    if !can_contract(i, j)
      msg = "Indices $i and $j cannot be contracted, see [`can_contract`](@ref). Argument sizes: $Sa and $Sb."
      return :( throw(DimensionMismatch($(msg) )) )
    end
  end

  S_contract = tuple(last_N_Sa...)
  Sa_keep, Sb_keep = sz_a[1:Na-N], sz_b[N+1:Nb]

  Sr = tuple(Sa_keep..., Sb_keep...)
  Nr = length(Sr)
  V = _tensor_base_type_expr(Nr,Sr)

  (iszero(length(a)) || iszero(length(b))) && return _zero_tensor_expr(V,Nr,Ta,Tb)

  # TODO, if the length of the resulting tensor exceeds some threshold, switch
  # to runtime looping over the indices, or even using BLAS. Also, we might need
  # to change NTensorValue to store into Memory or simply some fixed sized array.
  # Context: compiling double_contaction of 4th order 4D tensors takes 5-10 min, runtime 15μs.
  comps = Expr[]
  for cib in CartesianIndices(Int.(Sb_keep)) # Julia is column major, last index enumerates first
    for cia in CartesianIndices(Int.(Sa_keep))
      summands = [:(a[$cia, $ciC]*b[$ciC, $cib]) for ciC in CartesianIndices(Int.(S_contract))]
      push!(comps, Expr(:call, :+, summands...))
    end
  end
  iszero(Nr) ? :( $(comps[1]) ) : :( $V($(comps...)) )
end

# This method is for the material hodge stars: the most general way to transform
# a K form into a D-K form is to contract it with a (K-form)⊗(K-vector) tensor
@inline function contracted_product(
  ::Val{K},
  a::ProductTensorValue{A,B},
  b::ExteriorFormValue{K,Db}
  ) where {K,A,Da,B<:ExteriorVectorValue{K,Da},Db}


  Da == Db || throw(DimensionMismatch(
    "A $K-vector block of dimension $Da does not contract a $K-form of dimension $Db."))

  T = Base.promote_op(*, eltype(a), eltype(b))
  iszero(num_indep_components(A)) && return zero(make_concretetype(A{T}))
  make_concretetype(A{T})(_map_second_factor(x -> _flat_inner_body(x, b), a).data)
end

###############################################################
# General tensor contraction (arbitrary index pairs)
###############################################################

"""
    tensor_contraction(a::MultiValue, b::MultiValue, ia::Int, ib::Int)
    tensor_contraction(a::MultiValue, b::MultiValue, ia::NTuple{N,Int}, ib::NTuple{N,Int})
    tensor_contraction(a::MultiValue, b::MultiValue, ::Val{ia}, ::Val{ib})

Contract `N` index pairs between tensors `a` and `b`: the `ia[k]`-th index of `a`
is summed against the `ib[k]`-th index of `b`, for `k = 1, …, N`.
The contracted dimensions must match pairwise.

The output tensor has the remaining (non-contracted) indices of `a` in their
original order, followed by those of `b` in their original order.

!!! warning
    The methods accepting `Int`/`Tuple` are not type inferable and they allocate,
    pass the arguments by `Val` if they are known at compile time.

Generalises [`contracted_product`](@ref): `contracted_product(Val(N), a, b)` is
equivalent to `tensor_contraction(a, b, (Na-N+1, …, Na), (1, …, N))`.
"""
function tensor_contraction(
  a::MultiValue{Sa,Ta,Na},
  b::MultiValue{Sb,Tb,Nb},
  ia::NTuple{Nc,Int},
  ib::NTuple{Nc,Int}) where {Sa,Ta,Na,Sb,Tb,Nb,Nc}
  tensor_contraction(a, b, Val(ia), Val(ib))
end

tensor_contraction(a::MultiValue, b::MultiValue, ia::Int, ib::Int) =
  tensor_contraction(a, b, (ia,), (ib,))

@generated function tensor_contraction(
  a::MultiValue{Sa,Ta,Na},
  b::MultiValue{Sb,Tb,Nb},
  ::Val{Ia},
  ::Val{Ib}) where {Sa,Ta,Na,Sb,Tb,Nb,Ia,Ib}

  @assert (Ia isa NTuple || Ia isa Dim)
  @assert (Ib isa NTuple || Ib isa Dim)
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

  sz_a, sz_b = size(a), size(b)
  for k in 1:Nc
    da, db = sz_a[Ia[k]], sz_b[Ib[k]]
    if !can_contract(da, db)
      msg = "Index pair $k: Sa[$(Ia[k])]=$da and Sb[$(Ib[k])]=$db, cannot be contracted. See [`can_contract`](@ref)."
      return :(throw(DimensionMismatch($msg)))
    end
  end

  Ia_set = Set(Ia)
  Ib_set = Set(Ib)
  ka = tuple(filter(i -> !(i in Ia_set), 1:Na)...)
  kb = tuple(filter(i -> !(i in Ib_set), 1:Nb)...)

  S_a_kept   = ntuple(j -> sz_a[ka[j]], length(ka))
  S_b_kept   = ntuple(j -> sz_b[kb[j]], length(kb))
  S_contract = ntuple(k -> sz_a[Ia[k]], Nc)

  Sr = (S_a_kept..., S_b_kept...)
  Nr = length(Sr)
  V = _tensor_base_type_expr(Nr,Sr)

  (iszero(length(a)) || iszero(length(b))) && return _zero_tensor_expr(V,Nr,Ta,Tb)

  s_a_ranges = ntuple(j -> 1:Int(S_a_kept[j]), length(ka))
  s_b_ranges = ntuple(j -> 1:Int(S_b_kept[j]), length(kb))
  s_c_ranges = ntuple(k -> 1:Int(S_contract[k]), Nc)

  comps = Expr[]
  for cib in Iterators.product(s_b_ranges...)
    for cia in Iterators.product(s_a_ranges...)
      summands = Expr[]
      for ciC in Iterators.product(s_c_ranges...)
        a_idx = zero(MVector{Na,Int})
        for (j, pos) in enumerate(ka)
          a_idx[pos] = cia[j]
        end
        for (k, pos) in enumerate(Ia)
          a_idx[pos] = ciC[k]
        end
        b_idx = zero(MVector{Nb,Int})
        for (j, pos) in enumerate(kb)
          b_idx[pos] = cib[j]
        end
        for (k, pos) in enumerate(Ib)
          b_idx[pos] = ciC[k]
        end
        push!(summands, :(a[$(a_idx...)]*b[$(b_idx...)]) )
      end
      push!(comps, Expr(:call, :+, summands...))
    end
  end
  iszero(Nr) ? :( $(comps[1]) ) : :( $V($(comps...)) )
end

"""
    tensor_contraction(a::MultiValue, i::Int, j::Int)
    tensor_contraction(a::MultiValue, i::NTuple{N,Int}, j::NTuple{N,Int})
    tensor_contraction(a::MultiValue, ::Val{i}, ::Val{j})

Self-contraction of `a`: for each pair `k`, index `i[k]` and index `j[k]` of `a`
are set equal and summed over. `i` and `j` must be disjoint and the paired
dimensions must match.

The output has the remaining indices of `a` (those not in `i` or `j`) in their
original order. `tr(a)` is the special case `tensor_contraction(a, (1,), (2,))`.

!!! warning
    The methods accepting `Int`/`Tuple` are not type inferable and they allocate,
    pass the arguments by `Val` if they are known at compile time.
"""
function tensor_contraction(
  a::MultiValue{Sa,Ta,Na},
  i::NTuple{Nc,Int},
  j::NTuple{Nc,Int}) where {Sa,Ta,Na,Nc}
  tensor_contraction(a, Val(i), Val(j))
end

tensor_contraction(a::MultiValue, i::Int, j::Int) =
  tensor_contraction(a, (i,), (j,))

@generated function tensor_contraction(
  a::MultiValue{Sa,Ta,Na},
  ::Val{I},
  ::Val{J}) where {Sa,Ta,Na,I,J}

  @assert (I isa NTuple || I isa Int)
  @assert (J isa NTuple || J isa Int)
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

  sz_a = size(a)
  for k in 1:Nc
    di, dj = sz_a[I[k]], sz_a[J[k]]
    if !can_contract(di, dj)
      msg = "Index pair $k: Sa[$(I[k])]=$di and Sb[$(J[k])]=$dj, cannot be contracted. See [`can_contract`](@ref)."
      return :(throw(DimensionMismatch($msg)))
    end
  end

  contracted = Set([I..., J...])
  ka         = tuple(filter(idx -> !(idx in contracted), 1:Na)...)
  Sr     = ntuple(m -> sz_a[ka[m]], length(ka))
  S_contract = ntuple(k -> sz_a[I[k]], Nc)
  Nr         = length(Sr)
  V = _tensor_base_type_expr(Nr,Sr)

  iszero(length(a)) && return _zero_tensor_expr(V,Nr,Ta,Ta)

  s_keep_ranges     = ntuple(m -> 1:Int(Sr[m]), length(ka))
  s_contract_ranges = ntuple(k -> 1:Int(S_contract[k]), Nc)

  comps = Expr[]
  for cia in Iterators.product(s_keep_ranges...)
    summands = Expr[]
    for ciC in Iterators.product(s_contract_ranges...)
      a_idx = zero(MVector{Na,Int})
      for (m, pos) in enumerate(ka)
        a_idx[pos] = cia[m]
      end
      for (k, pos) in enumerate(I)
        a_idx[pos] = ciC[k]
      end
      for (k, pos) in enumerate(J)
        a_idx[pos] = ciC[k]
      end
      push!(summands, :(a[$(a_idx...)]))
    end
    push!(comps, Expr(:call, :+, summands...))
  end
  iszero(Nr) ? :( $(comps[1]) ) : :( $V($(comps...)) )
end

###############################################################
# Congruent product
#
# bᵀ⋅a⋅b if a is order 2,
#
# a_J = Σ_I M[j₁,i₁] ⋯ M[j_N,i_N] a_I else.
#
###############################################################

# tensor type of the result of congruent_prod preserving as much symmetries as
# possible, depending on the argument types
_congruence_type(::Type{Va}, ::Type{Vb}) where {Va<:MultiValue,Vb<:MultiValue} =
  _congruence_type(Va, size(Vb)[2])

_congruence_type(::Type{V}, Dout::Dim) where V<:MultiValue =
  _tensor_base_type(length(size(V)), ntuple(_ -> Dout, length(size(V))))

_congruence_type(::Type{<:AbstractSymTensorValue}, Dout::Dim)    = SymTensorValue{Dout}
_congruence_type(::Type{<:SkewSymTensorValue}, Dout::Dim)        = SkewSymTensorValue{Dout}
_congruence_type(::Type{<:SymFourthOrderTensorValue}, Dout::Dim) = SymFourthOrderTensorValue{Dout}
_congruence_type(::Type{V}, Dout::Dim) where V<:SymNTensorValue =
  SymNTensorValue{Tuple{ntuple(_ -> Dout, length(size(V)))...}}
_congruence_type(::Type{<:SkewSymKTensorValue{V,K}}, Dout::Dim) where {V,K} =
  SkewSymKTensorValue{typeof(Dout),K,Dout}

"""
    congruent_prod(a, b)
    congruent_prod(a, b, ::Val{l})

Given an order `N ≥ 1` tensor `a` and a second order tensor `b`, return the
congruence of every index of `a` by `b`

    c_O = Σ_I b[i₁,o₁] ⋯ b[i_N,o_N] a_I

that is `b`ᵀ⋅`a`⋅`b` for `N = 2`, and `b`ᵀ⋅`a` for `N = 1`.
The type of the resulting value is (skew) symmetric stable w.r.t. `typeof(a)`,
a traceless `a` widening to a `SymTensorValue`: a congruence keeps an index
symmetry but not a vanishing trace, see [`change_variance`](@ref) for the
relabelling that keeps both.

`a` carries `N` times the same index `D`, variance included, and `b`ᵀ contracts
its second index against each of them, so the first index `Db` of `b` must
satisfy `can_contract(Db,D)`. Every index of the result carries the second
index of `b`, and only the components the result stores are computed.

A congruence needs an index to contract, so an order 0 `a` is an error.

The three argument form congruences the indices `k ∈ l` only,

    c_{…,o_k,…} = Σ_{…,i_k,…} b[i_k,o_k] ⋯ a_{…,i_k,…}

leaving the others untouched. Only the indices in `l` need to share their
dimension and variance, the others are unconstrained. An index symmetry
survives a congruence of *every* index, so a `l` that leaves one out returns a
tensor of the corresponding full rank type: the widening is automatic, it is
the two argument form that preserves symmetries.
"""
function congruent_prod(a::Number, b::Number)
  msg = """ operation only defined for 2nd order tensors `a` and `b` with
      `size(b,1) == size(a, 1) == size(a, 2)`, got `size(a)=$(size(a))` and
      `size(b)=$(size(b))`.
      """
  @unreachable msg
end

congruent_prod(a::Number, b::Number, ::Val) = congruent_prod(a, b)

# Contract the index `k` of tensor data `x`, of dimension `Din`, into the index
# `k` of `y`, of dimension `Dout`
#
#   y[i1,…,ik,…,iN] = Σ_j M[j,ik] * x[i1,…,j,…,iN]
#
# `M`'s contracted index comes first and its free one gives the result, the
# orientation of the `b` of a congruence `b`ᵀ⋅`a`⋅`b`. A metric is symmetric, so
# this orientation also lowers or raises an index.
#
# `y` and `x` must be distinct buffers. The stride `s` of `k` is the number of
# components between two that differ only by ±1 their `k`th index, and is shared
# by the two buffers: the indices before `k` have the same dimension in both, only
# `k` and those after it may differ. The indices after `k` give the
# `Lin ÷ (s*Din) == Lout ÷ (s*Dout)` blocks the contraction runs over.
@inline function _data_single_congruence_contract!(
  y::MVector{Lout,T}, x::StaticVector{Lin}, M, ::Val{s}, ::Val{Din}, ::Val{Dout}
) where {Lout,Lin,T,s,Din,Dout}

  # TODO: check in marking the whole block @inbounds or @simd makes it faster
  for o in 0:(Lin ÷ (s*Din) - 1)
    for p in 1:s
      qin, qout = o*s*Din + p, o*s*Dout + p
      for i in 1:Dout
        acc = zero(T)
        for j in 1:Din
          acc += M[j,i] * x[qin + (j-1)*s]
        end
        y[qout + (i-1)*s] = acc
      end
    end
  end
  nothing
end

# number of tensor terms above which the congruence operations are not all
# hard-coded in a single @generated method body ("direct"), but split in
# several steps.
const _congruence_direct_term_threshold = 200

# Index "O with ks ← J", that is O but with O[i] ← J[m] if i = ks[m]
_congruence_source_index(O, ks, J) =
  ntuple(i -> (m = findfirst(==(i), ks)) === nothing ? O[i] : J[m], length(O))

# `c_O = Σ_J a_{O with ks ← J} ∏ₘ b[J[m], O[ks[m]]]`, the fully unrolled
# congruence over the indices `ks`, of dimension `Din`
function _congruence_direct_body(Vr, T, ks, Din)
  Js = CartesianIndices(ntuple(_ -> Din, length(ks)))
  comps = [ _sum_expr(
              [ Expr(:call, :*, :( a[$(_congruence_source_index(O, ks, J)...)] ),
                     (:( b[$(J[m]),$(O[ks[m]])] ) for m in eachindex(ks))...)
                for J in Js ], T)
            for O in _stored_comp_cartesian_indices(make_concretetype(Vr{Int})) ]
  :( $Vr(($(comps...),)) )
end

# `c_O = Σ_I b[i₁,o₁] … b[i_N,o_N] a_I`, the congruence of every index of `a` by
# the one matrix `b`, which has as many rows as `a` has components along one
# index.
#
# The `N-1` first indices are contracted by the mode-product ping-pong the
# generic methods use, whose loop LLVM vectorizes. The last one is contracted
# straight into the components the result stores, which costs
# `num_indep_components·Din` products where one more iteration would cost
# `Dout^N·Din`.
#
# Only reading the stored components was tried and didn't execute faster.
@generated function congruent_prod(a::MultiValue{S,Ta,N},
                                   b::MultiValue{Tuple{Db,D1},Tb}) where {S,Ta,N,Db,D1,Tb}

  sa = size(a)

  msg = """ operation only defined for a tensor `a` of order at least 1 whose
      indices are all identical, variance included, and a 2nd order `b` with
      `can_contract(size(b,1), size(a,1))`, got `size(a)=$sa` and
      `size(b)=$((Db,D1))`.
      """
  # Order 0 carries no index, so it is refused before `sa` is read.
  iszero(N) && return :( @unreachable $msg )

  Da = first(sa)
  all(D -> D === Da, sa) || return :( @unreachable $msg )
  can_contract(Db, Da) || return :( _ensure_contract($Db, $Da) )

  Vr = _congruence_type(a, b)
  Din, Dout = Int(Da), Int(D1)
  T = :( Base.promote_op(*, $Tb, $Ta) )
  (iszero(Din) || iszero(Dout)) && return :( zero($Vr{$T}) )

  ks = ntuple(identity, N)
  num_terms = num_stored_components(make_concretetype(Vr{Int})) * Din^N
  num_terms <= _congruence_direct_term_threshold && return _congruence_direct_body(Vr, T, ks, Din)

  buf_len(k) = Dout^(k-1) * Din^(N-k+1)
  if Din == Dout
    bufs = ntuple(k -> isodd(k) ? :x0 : :x1, N)
    define_buffers = N > 1 ? Expr(:block, :( x1 = MVector{$(buf_len(1)),$T}(undef) )) : Expr(:block)
  else
    bufs = ntuple(k -> Symbol(:x, k-1), N)
    define_buffers = Expr(:block,
                 (:( $(bufs[k]) = MVector{$(buf_len(k)),$T}(undef) ) for k in 2:N)...)
  end
  congruence_contractions = Expr(:block,
    (:( _data_single_congruence_contract!($(bufs[k+1]), $(bufs[k]), b,
          Val($(Dout^(k-1))), Val($Din), Val($Dout)) ) for k in 1:N-1)...)

  # Position of the component `(O₁, …, O_{N-1}, j)` of the buffer the
  # congruence_contractions leave the `N-1` first indices converted in, the
  # last one having stride `Dout^(N-1)`.
  pos(O, j) = 1 + sum(((O[k]-1)*Dout^(k-1) for k in 1:N-1), init=0) + (j-1)*Dout^(N-1)
  comps = [ _sum_expr([ :( b[$j,$(O[N])] * y[$(pos(O,j))] ) for j in 1:Din ], T)
            for O in _stored_comp_cartesian_indices(make_concretetype(Vr{Int})) ]

  quote
    x0 = MVector{$(buf_len(1)),$T}(Tuple(get_array(a)))
    $define_buffers
    $congruence_contractions
    y = Tuple($(bufs[N])) # for more performent getindex in comps access
    $Vr(($(comps...),))
  end
end

# The congruence_contractions of a congruence over the indices `ks` of a tensor
# of size `sz`, each contracted from its dimension `Dk=sz[k]` to `Dout`.
# Returns for every `k ∈ `ks`:
# - stride: the stride that `k` has during it's congruence contraction,
# - len_in: length of the input buffer of that iteration,
# - len_out:length of the output buffer ..
#
# The length change at iteration k happens if Dout ≠ Dk
function _congruence_contractions_data(sz, ks, Dout)
  dims = [Int(d) for d in sz]
  map(ks) do k
    s, len_in = prod(dims[1:k-1], init=1), prod(dims, init=1)
    dims[k] = Dout
    (stride = s, len_in = len_in, len_out = prod(dims, init=1))
  end
end

# `c_{…o_k…} = Σ_{…i_k…} b[i_k,o_k] ⋯ a_{…i_k…}` over the indices `k ∈ l` only.
# preserving symmetries if possible
@generated function congruent_prod(a::MultiValue{S,Ta,N},
    b::MultiValue{Tuple{Db,D1},Tb}, ::Val{l}) where {S,Ta,N,Db,D1,Tb,l}

  ks = tuple(l...) # puts a scalar `l` into a tuple
  isempty(ks) && return :( a )

  sz = size(a)
  # Order 0 carries no index, so it is refused before `sz` is read.
  iszero(N) && return :( @unreachable $(
    "congruent_prod is only defined for a tensor `a` of order at least 1, got order 0") )

  all(k -> k isa Int && 1 <= k <= N, ks) || return :( throw(ArgumentError(
    $("congruent_prod undefined for the indices $ks of a tensor of order $N"))) )
  allunique(ks) || return :( throw(ArgumentError(
    $("congruent_prod cannot act twice on one index, got the indices $ks"))) )

  # Every index congruenced: the symmetry preserving body applies.
  length(ks) == N && return :( congruent_prod(a, b) )

  Da = sz[first(ks)]
  all(k -> sz[k] === Da, ks) || return :( throw(DimensionMismatch(
    $("the indices $ks of a tensor of size $sz are not all of the same dimension \
       and variance, so no single `b` can congruence them all"))) )
  can_contract(Db, Da) || return :( _ensure_contract($Db, $Da) )

  Sr = ntuple(i -> i in ks ? D1 : sz[i], N)
  Vr = _tensor_base_type(N, Sr)
  T = :( Base.promote_op(*, $Tb, $Ta) )
  Din, Dout = Int(Da), Int(D1)

  # zero output
  (iszero(length(a)) || iszero(prod(Int(d) for d in Sr; init=1))) && return :( zero($Vr{$T}) )

  # if small number of operations, unroll them all in one body
  num_terms = num_stored_components(make_concretetype(Vr{Int})) * Din^length(ks)
  num_terms <= _congruence_direct_term_threshold && return _congruence_direct_body(Vr, T, ks, Din)

  # else, compute the contraction for each index in a separately compiled body
  ps = _congruence_contractions_data(sz, sort!(collect(ks)), Dout)
  nk = length(ps)
  lens = (first(ps).len_in, (p.len_out for p in ps)...)
  if Din == Dout
    bufs = ntuple(m -> isodd(m) ? :x0 : :x1, nk+1)
    define_buffers = Expr(:block, :( x1 = MVector{$(lens[1]),$T}(undef) ))
  else
    bufs = ntuple(m -> Symbol(:x, m-1), nk+1)
    define_buffers = Expr(:block,
      (:( $(bufs[m]) = MVector{$(lens[m]),$T}(undef) ) for m in 2:nk+1)...)
  end
  congruence_contractions = Expr(:block,
    (:( _data_single_congruence_contract!($(bufs[m+1]), $(bufs[m]), b,
          Val($(ps[m].stride)), Val($Din), Val($Dout)) ) for m in 1:nk)...)

  quote
    x0 = MVector{$(lens[1]),$T}(Tuple(get_array(a)))
    $define_buffers
    $congruence_contractions
    $Vr(Tuple($(bufs[nk+1])))
  end
end

# `c_O = Σ_{I sorted} det(b[I,O]) a_I`, the congruence of every index of an
# alternating tensor of order `K` by `b`, from dimension `Din` to `Dout`.
function _alt_congruence_body(Vr, K, Din, Dout, T)
  (iszero(binomial(Din, K)) || iszero(binomial(Dout, K))) && return :( zero($Vr{$T}) )
  comps = _alt_compound_comps(:b, sorted_combinations(Dout, K),
                              sorted_combinations(Din, K), T, :a, true)
  :( $Vr(($(comps...),)) )
end

@generated function congruent_prod(a::SkewSymKTensorValue{V,K,Dw,Ta},
    b::MultiValue{Tuple{Db,D1},Tb}) where {V,K,Dw,Ta,Db,D1,Tb}

  can_contract(Db, Dw) || return :( _ensure_contract($Db, $Dw) )
  T = Base.promote_op(*, Tb, Ta)
  Vr = _congruence_type(a, b)
  iszero(K) && return :( $Vr{$T}(a.data) )
  _alt_congruence_body(Vr, K, Int(Dw), Int(D1), T)
end

# same as `SkewSymKTensorValue{V,2,D}`
@generated function congruent_prod(a::SkewSymTensorValue{Dw,Ta},
    b::MultiValue{Tuple{Db,D1},Tb}) where {Dw,Ta,Db,D1,Tb}

  can_contract(Db, Dw) || return :( _ensure_contract($Db, $Dw) )
  _alt_congruence_body(_congruence_type(a, b), 2, Int(Dw), Int(D1),
                       Base.promote_op(*, Tb, Ta))
end
