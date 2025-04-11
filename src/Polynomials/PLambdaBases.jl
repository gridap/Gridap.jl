
#####################
# Combination type  #
#####################

import Base

const MAX_PRECOMP_D = 8
const NB_PRECOMP_IPERMS = 2^MAX_PRECOMP_D
const Val_Precomp_D = Union{Val{0},Val{1},Val{2},Val{3},Val{4},Val{5},Val{6},Val{7},Val{8}}

"""
Combination{k,D} <: AbstractSet{Int}


"""
struct Combination{k,D} <: AbstractSet{Int}
  data::NTuple{k,Int}

  function Combination{k,D}(perm::NTuple{k,Int}) where {k,D}
    @check 0 ≤ k ≤ D "0 ≤ k ≤ D is required, got k=$k and D=$D"
    @check issorted(perm) && allunique(perm) "the given tuple is not (strictly) increasing, perm=$perm"
    @check all(@. 1 ≤ perm ≤ D) "the given tuple is not a subset of ⟦1,`D`⟧, perm=$perm"
    new{k,D}(perm)
  end
end

Base.convert(::Type{NTuple{k,Int}}, x::Combination{k}) where k = x.data

Combination(perm, D) = Combination(perm, Val(D))
Combination(perm, ::Val{D}) where D = Combination{length(perm),D}(Tuple(sort(unique(perm))))
Combination{0,D}() where D = Combination{0,D}( () )
function Combination{k,D}(v::AbstractVector) where {k,D}
  @check k == length(v) "Isnconsistentcy between permutation length and number of given indices"
  Combination{k,D}( Tuple(v) )
end

"""
_ID(combi::Combination{k})
_ID(data::NTtuple{k,Int})

Internal iddentifyer of a permutation of length `k` (it is independant of the second struct parameter).
"""
_ID(combi::Combination) = _ID(combi.data)
_ID(data::NTuple{k,Int}) where k = mapreduce(t -> 1<<(t-1), +, data, init=0) + 1

# Not type stable, do not use at runtime
function _ID_to_data(ID::Integer)
  @check ID > 0
  bits = ID-1
  curr, pos = 0, 1
  k = count_ones(bits)
  data = MVector{k,Int}(undef)
  while curr < k
    if bits & 1 ≠ 0
      curr += 1
      data[curr] = pos
    end
    bits = bits >> 1
    pos += 1
  end
  Tuple(data)
end

_ID_to_k(ID::Integer) = count_ones(ID-1)


# Iteration Interface
Base.length(::Combination{k}) where k = k
Base.iterate(::Combination{0}) = nothing
Base.isdone(combi::Combination) = Base.isdone(combi.data)
Base.isdone(combi::Combination, state) = Base.isdone(combi.data, state)

Base.iterate(combi::Combination) = iterate(combi.data)
Base.iterate(combi::Combination, state) = iterate(combi.data, state)

# Indexing Interface
Base.IndexStyle(::Combination) = IndexLinear()
Base.keys(::Combination{k}) where k = Base.OneTo(k)
Base.firstindex(::Combination) = 1
Base.lastindex(::Combination{k}) where k = k

function Base.getindex(combi::Combination{k}, i) where k
  @check 1 ≤ i ≤ k "Indices in Combinations{$k} are 1:$k, got $i"
  combi.data[i]
end

# Iterable Collection API

Base.in(item, ipa::Combination) = in(item, ipa.data)
Base.indexin(items, combi::Combination) = indexin(items, combi.data)

Base.maximum(combi::Combination{k}; init) where k = combi[k]
Base.maximum(combi::Combination{0}; init) = maximum((); init=Int(init))

Base.minimum(combi::Combination{k}; init) where k = combi[1]
Base.minimum(combi::Combination{0}; init) = minimum((); init=Int(init))


# Set & Ordering API

for fun in (:union, :intersect, :setdiff, :symdiff, :issubset, :issetequal, :isdisjoint, :isequal)
  @eval begin
    Base.$fun(a::Combination, b::Combination) = $fun(a.data, b.data)
  end
end

Base.isunordered(combi::Combination) = false

"""
Base.isless(a::Combination, b::Combination)

The ordering of `Combination` is right-digit to left-digit
lexicographic order, e.g. 12 < 13 < 23 < 14 < 24 < 34.

Unlike with the usual (left-digit to right-digit) lexicographic order, the
index of sorted `Combination`s of same length `k` is independant of the
dimension `D`, see [`sorted_combinations`](@ref).
"""
function Base.isless(a::Combination, b::Combination)
  isless(reverse(a.data), reverse(b.data))
end

# Other API

"""
sorted_combinations(D,k)
sorted_combinations(::Val{D},::Val{k})

Return a vector [Iᵢ]ᵢ of all the combinations of 1:D of length k:

1≤ I_1 < ... < I_k ≤ N

sorted in right-to-left lexicographic order, e.g.

[12, 13, 23]               for k=2, D=3
[12, 13, 23, 14, 24, 34]   for k=2, D=4

See also [`Base.isless`](@ref).
"""
sorted_combinations(D,k) = sorted_combinations(Val(D),Val(k))

@generated function sorted_combinations(::Val{D},::Val{k}) where {D,k}
  iszero(k) && return :( return [Combination{0,D}( () )] )
  comp_rev_perm(tup) = ntuple(i -> D-tup[k-i+1]+1, Val(k))
  inc_perms = combinations(1:D,k) .|> (tup -> Combination{k,D}(comp_rev_perm(tup))) |> reverse
  :( return $inc_perms )
end

"""
complement(combi::Combination{k,D}) where {k,D}

TBW
"""
function complement(combi::Combination{k,D}) where {k,D}
  comp = MVector{D-k,Int}(undef)
  curr_perm, curr_comp = 1, 1
  for i in 1:D
    if curr_perm ≤ k && combi[curr_perm] == i
      curr_perm += 1
    else
      comp[curr_comp] = i
      curr_comp += 1
    end
  end
  return Combination{D-k,D}(Tuple(comp))
end

"""
combination_sign(combi::Combination)

TBW
"""
function combination_sign(combi::Combination)
  i, k, acc, delta = 1, 1, 0, 0
  while k <= length(combi)
    if combi[k] == i
      acc += delta
      k += 1
    else
      delta += 1
    end
    i += 1
  end
  return iseven(acc) ? 1 : -1
end

"""
combination_index(combi::Combination)

Linear index of `combi` amongst combinations of the same size `k`,
sorted in right-to-left lexicographic order. It depends on `k` but not on `D`,
see [`sorted_combinations`](@ref).
"""
function combination_index(combi::Combination{k,D}) where {k,D}
  combination_index(combi, Val(D))
end
@inline function combination_index(combi::Combination, ::Val_Precomp_D)
  precomp_ID_to_combi_index[_ID(combi)]
end

@inline function combination_index(combi::Combination{k,D}, _) where {k,D}
  ID = _ID(combi)
  if haskey(memo_ID_to_combi_index, ID)
    return memo_ID_to_combi_index[ID]
  else
    combi_index = findfirst(==(combi), sorted_combinations(Val(D),Val(k)))
    isnothing(combi_index) && @unreachable
    memo_ID_to_combi_index[ID] = combi_index
    return combi_index
  end
end

const memo_ID_to_combi_index = IdDict{Int,Int}()

function generate_ID_to_combi_index()
  ID_to_combi_index = MVector{NB_PRECOMP_IPERMS,Int}(undef)
  curr_indices = zero(MVector{MAX_PRECOMP_D+1,UInt8})
  for ID in 1:NB_PRECOMP_IPERMS
    k = _ID_to_k(ID)
    curr_indices[k+1] += 1
    ID_to_combi_index[ID] = curr_indices[k+1]
  end
  SVector(ID_to_combi_index)
end

const precomp_ID_to_combi_index = generate_ID_to_combi_index()


"""
sub_combinations(combi::Combination{k,D})

Return a tuple containing the `k-1` combinations `combi\\combi[i]` for 1 ≤ i ≤ `k`.
"""
sub_combinations(::Combination{0,D}) where D = tuple()
function sub_combinations(combi::Combination{k,D}) where {k,D}
  @check k>0
  iszero(k) && return Tuple{}[]
  m1_combi_datas = MVector{k,Combination{k-1,D}}(undef)
  for i in 1:k
    mi_data = ntuple(j -> combi.data[j + Int(j≥i)],Val(k-1))
    m1_combi_datas[i] = Combination{k-1,D}(mi_data)
  end
  return Tuple(m1_combi_datas)
end

# Display
function Base.show(io::IO, combi::Combination{k,D}) where {k,D}
  iszero(k) && print(io,"∅")
  print(io,"$(join(combi.data))")
end
function Base.show(io::IO, ::MIME"text/plain", combi::Combination{k,D}) where {k,D}
  if iszero(k)
    print(io,"Combination{$k,$D}(∅)")
  else
    print(io,"Combination{$k,$D}($(join(combi.data,",")))")
  end
end
#################################
# Tensorial nD polynomial bases #
#################################

"""
struct PLambda{D,T,L,K,B,P} <: PolynomialBasis{D,VectorValue{L,T},K,Bernstein}

Finite Element Exterior Calculus polynomial basis for the spaces `D`-dimensional
simplices, P`r`Λ`ᴷ`(△`ᴰ`).

Reference: D. N. Arnold and A. Logg, Periodic Table of the Finite Elements, SIAM News, vol. 47 no. 9, November 2014
"""
struct PLambdaBasis{D,T,L,K,P,B} <: PolynomialBasis{D,VectorValue{L,T},K,Bernstein}
  r::Int
  k::Int
  Ψ::P
  b::B

  function PLambdaBasis{D}(::Type{T},r,k,vertices=nothing) where {D,T}
    @check T<:Real "T needs to be <:Real since represents the scalar type"
    @check k in 0:D "The form order k must be in 0:D"
    @check r > 0    "The polynomial order r must be positive"
    if !isnothing(vertices)
      @check length(vertices) == D+1 "$D+1 vertices are required to define a $D-dim simplex, got $(length(vertices))"
      @check eltype(vertices) isa Point{D} "Vertices should be of type Point{$D}, got $(eltype(vertices))"
    end

    L = binomial(k,D) # Number of components of a basis form
    C = binomial(r+k,k)*binomial(D+r,D-k) # Number of basis polynomials

    b = BernsteinBasisOnSimplex{D}(T, r, vertices=vertices)
    B = typeof(b)

    if isnothing(vertices)
      Ψ = nothing
      P = Nothing
    else
      # Col major: the L components of a basis form are stored contiguously
      P = SMatrix{L,C,T}
      Ψ = MMatrix{L,C,T}(undef)
      Ψ = _compute_PΛ_basis_form_coefficient!(Ψ,r,k,D,b,vertices)
      # Ψ = convert(P, Ψ) # isn't this automatic ?
    end

    new{D,T,L,K,P,B}(r,k,Ψ,b)
  end
end

function PLambdaBasis(::Val{D},::Type{T},r,k,vertices=nothing) where {D,T}
  PLambdaBasis{D}(T,r,vertices)
end

get_FEEC_poly_degree(b::PLambdaBasis) = b.r
get_FEEC_form_degree(b::PLambdaBasis) = b.k
get_FEEC_family(::PLambdaBasis) = :P

function Base.size(b::PLambdaBasis{D}) where D
  r, k = b.r, b.k
  return (binomial(r+k,k)*binomial(D+r,D-k), )
end

Base.getindex(b::PLambdaBasis, i::Integer) = Bernstein()
Base.IndexStyle(::PLambdaBasis) = IndexLinear()
get_dimension(::PLambdaBasis{D}) where D = D
get_order(b::PLambdaBasis) = get_FEEC_poly_degree(b)
return_type(::PLambdaBasis{D,T,L}) where {D,T,L} = VectorValue{L,T}


###################
# Implementation  #
###################

_compute_PΛ_basis_form_coefficient!(Ψ,r,k,D,b,::Nothing) = nothing # done in _hat_Ψ

function _compute_PΛ_basis_form_coefficient!(Ψ,r,k,D,b,vertices)
  N = D+1
  T = eltype(Ψ)
  α_prec = ntuple(_->-1, N)
  φ_αF = MMatrix{D,N,T}(undef)
  for (_, F, dF_bubbles) in PΛ_bubble_indices(r,k,D)
    for (w, α, J) in dF_bubbles
      if α ≠ α_prec
        update_φ_αF!(φ_αF,b,α,F,r)
        α_prec = α
      end
      for (iI,I) in enumerate(sorted_combinations(Val{k},Val{D}))
        Ψ[iI,w] = @inline minor(φ_αF,I,J)
      end
    end
  end
end

@inline function update_φ_αF!(φ_αF,b,α,F,r)
  M = b.cart_to_bary_matrix
  @inbounds for ci in CartesianIndices(φ_αF)
    i, j = ci[1], ci[2]
    mF = sum(M[Fl,i+1] for Fl in F; init=0)
    φ_αF[i][j] = M[j,i+1] - α[j]*mF/r
  end
end

"""
_hat_Ψ(r,α,F,I,J,T)::T

PLambdaBasis.Ψ matrix elements in the reference simplex, T is the scalar return type
"""
function _hat_Ψ(r,α,F,I,J,::Type{T})::T where T
  @check sum(α) == r
  @check length(I) == length(J) # thanks to dispatch on zero length J below
  @check length(J) > 0 # thanks to dispatch on zero length J below

  @inbounds begin

    s = Int(isone(J[1]))
    n = count(i-> (J[i]-1)∉I, (s+1):length(J))

    n > 1 && return 0. # rank M_IJ inferior to 2

    p = _findfirst_or_zero(j-> (I[j]+1)∉J, 1:length(J))
    if iszero(p)
      println("p0 $I $J $s $n")
    end

    if isone(n)        # rank M_IJ is 1
      m = _findfirst_or_zero(i-> (J[i]-1)∉I, (s+1):length(J))
      if iszero(m)
        println("m0 $I $J $s $n")
      end
      u_m, v_p = _u(m,F,I), _v(p,α,J,r)
      iszero(s) && return (-1)^(m+p)*u_m*v_p

      q = _findfirst_or_zero(j-> (I[j]+s)∉J, (p+1):length(J))
      @check !iszero(q)
      v_q = _v(q,α,J,r)
      return (-1)^(m+p+q) * u_m * (v_q - v_p)
    end

    u, v = _u(F,I), _v(α,J,r)
    if iszero(s)
      return 1 + sum( u .* v )
    else
      Ψ_IJ = one(T)
      sum_v = v[1]
      for l in 1:p-1
        vlp = v[l+1]
        sum_v += vlp
        Ψ_IJ += vlp*u[l]
      end
      for l in (p+1):length(J)
        vl = v[l]
        sum_v += vl
        Ψ_IJ += vl*u[l]
      end
      return (-1)^p * (Ψ_IJ - u[p]*sum_v)
    end

  end
  @unreachable
end
_hat_Ψ(r,α,F,I,J::Combination{0},::Type{T}) where T = one(T) # 0 forms

function _findfirst_or_zero(pred, v)
  r = findfirst(pred,v)
  return isnothing(r) ? 0 : r
end

@propagate_inbounds _u(i,F,I)   = Int(isone(F[1])) - Int(I[i]+1 in F)
@propagate_inbounds _u(F,I::Combination{k}) where k = ntuple(i->_u(i,F,I), Val(k))
@propagate_inbounds _v(j,α,J,r) = α[J[j]]/r
@propagate_inbounds _v(α,J::Combination{k},r) where k = ntuple(j->_v(j,α,J,r), Val(k))


# API

function _return_cache(
  b::PLambdaBasis{D},x,::Type{G},::Val{N_deriv}) where {D,G,N_deriv}

  @assert D == length(eltype(x)) "Incorrect number of point components"
  T = eltype(G)
  r = get_order(b)
  np = length(x)
  ndof = length(b)
  ndof_bernstein = _binomial(Val(r+D),Val(D))
  ndof_val = _binomial(Val(r+D),Val(D))

  r = CachedArray(zeros(G,(np,ndof)))
  s = MArray{Tuple{Vararg{D,N_deriv}},T}(undef)
  # Cache for all scalar nD-Bernstein polynomials, no other caches needed for derivatives
  c = CachedVector(zeros(T,ndof_scalar))
  # Cache for 1
  t = ntuple( _ -> nothing, Val(N_deriv))
  (r, s, c, t...)
end

#function return_cache(b::PLambdaBasis, x::AbstractVector{<:Point})
#  return_cache(b._basis,x)
#end

#function evaluate!(cache, b::PLambdaBasis, x::AbstractVector{<:Point})
#  evaluate!(cache, b._basis, x)
#end

# Derivatives, option 2:

function return_cache(
  fg::FieldGradientArray{1,<:PLambdaBasis},
  x::AbstractVector{<:Point})

  fg_b = FieldGradientArray{1}(fg.fa._basis)
  return (fg_b, return_cache(fg_b,x))
end

function evaluate!(cache,
                   fg::FieldGradientArray{1,<:PLambdaBasis},
                   x::AbstractVector{<:Point})

  fg_b, b_cache = cache
  evaluate!(b_cache, fg_b, x)
end

function return_cache(
  fg::FieldGradientArray{2,<:PLambdaBasis},
  x::AbstractVector{<:Point})

  fg_b = FieldGradientArray{2}(fg.fa._basis)
  return (fg_b, return_cache(fg_b,x))
end

function evaluate!(cache,
                   fg::FieldGradientArray{2,<:PLambdaBasis},
                   x::AbstractVector{<:Point})

  fg_b, b_cache = cache
  evaluate!(b_cache, fg_b, x)
end


##########################
# PLambda bases helpers  #
##########################


PΛ_bubble_indices(r,k,D) = PΛ_bubble_indices(Val(r),Val(k),Val(D))
@generated function PΛ_bubble_indices(::Val{r},::Val{k},::Val{D}) where {r,k,D}
  i=0
  d_F_bubbles = []
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      dF_bubbles = _PΛ_F_bubble_indices(r,k,D,F,i)
      push!(d_F_bubbles, (d, F, dF_bubbles))
      i += length(dF_bubbles)
    end
  end
  @check i == binomial(r+k,k)*binomial(D+r,D-k)
  d_F_bubbles = tuple(d_F_bubbles...)
  :( $(d_F_bubbles) )
end

P⁻Λ_bubble_indices(r,k,D) = P⁻Λ_bubble_indices(Val(r),Val(k),Val(D))
@generated function P⁻Λ_bubble_indices(::Val{r},::Val{k},::Val{D}) where {r,k,D}
  i=0
  d_F_bubbles = []
  for d in k:D
    for F in sorted_combinations(D+1, d+1)
      dF_bubbles = _P⁻Λ_F_bubble_indices(r,k,D,F,i)
      push!(d_F_bubbles, (d, F, dF_bubbles))
      i += length(dF_bubbles)
    end
  end
  @check i == binomial(r+k-1,k)*binomial(D+r,D-k)
  d_F_bubbles = tuple(d_F_bubbles...)
  :( $(d_F_bubbles) )
end

function _P⁻Λ_F_bubble_indices(r,k,D,F,i)
  N = D + 1
  ids = []
  for α in bernstein_terms(r-1,N)
    for J in sorted_combinations(N,k+1)
      j = minimum(J)-1
      if issetequal(supp(α) ∪ J, F) && all(α[1:j] .== 0)
        push!(ids, (i, α, J))
        i += 1
      end
    end
  end
  tuple(ids...)
end

function _PΛ_F_bubble_indices(r,k,D,F,i)
  N = D + 1
  ids = []
  for α in bernstein_terms(r,N)
    for J in sorted_combinations(N,k)
      j = minimum(setdiff(F,J))-1
      if issetequal(supp(α) ∪ J, F) && all(α[1:j] .== 0)
        push!(ids, (i, α, J))
        i += 1
      end
    end
  end
  tuple(ids...)
end


"""
supp(α)

TBW
"""
function supp(α)
  s = Int[]
  for (i,αi) in enumerate(α)
    if αi > 0
      push!(s, i)
    end
  end
  Tuple(s)
end

function minor(M,I,J)
  @check length(I) == length(J)
  @check I ⊆ axes(M)[1]
  @check J ⊆ axes(M)[2]

  k = length(I)
  T = eltype(M)
  m = MMatrix{k,k,T}(undef)
  for (i, Ii) in enumerate(I)
    for (j, Jj) in enumerate(J)
      @inbounds m[i,j] = M[Ii,Jj]
    end
  end
  det(m)
end

function all_k_minors!(m,M,::Val{k}) where {k}
  D = size(M)[1]
  Λᵏᴰ = sorted_combinations(Val(D),Val(k))
  @inbounds begin
    for (i, I) in enumerate(Λᵏᴰ)
      for (j, J) in enumerate(Λᵏᴰ)
        m[i,j] = @inline minor(M,I,J)
      end
    end
  end
end

