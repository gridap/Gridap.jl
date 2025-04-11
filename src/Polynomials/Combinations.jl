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
Base.maximum(combi::Combination{0}; init) = Int(maximum((); init=init))

Base.minimum(combi::Combination{k}; init) where k = combi[1]
Base.minimum(combi::Combination{0}; init) = Int(minimum((); init=init))


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
