#####################
# Combination type  #
#####################

import Base

const Combination = Vector{Int}
const BubbleIndexSet = Tuple{Int, BernsteinTerm, Int, Combination}
const Bubble = Tuple{Int, Combination, Vector{BubbleIndexSet}}

struct Bubbles <: AbstractVector{Bubble}
  bubbles::Vector{Bubble}
end
Base.size(b_set::Bubbles) = size(b_set.bubbles)
Base.IndexStyle(::Type{Bubbles}) = IndexStyle(Vector{Bubble})
Base.getindex( b_set::Bubbles,      i::Int) = getindex( b_set.bubbles,      i::Int)
Base.setindex!(b_set::Bubbles, val, i::Int) = setindex!(b_set.bubbles, val, i::Int)

struct PLambdaIndices
  # An objectid can be stored here to keep track of what's inside and perform
  # sanity check when reusing the indices elsewhere, e.g. objectid( (r,k,D,:P) )
  # for a PᵣΛᵏ(△ᴰ) basis indices.
  identity::UInt

  # bubble indices given by P(m)Λ_bubbles
  bubbles::Bubbles

  # combination_index(I) -> (combination_index(I), I, sgnIcomp) ∀ I ∈ sorted_combinations(D,k)
  components::Vector{Tuple{Int, Combination, Int}}

  # _ID(J) -> ( l, _ID( J\{J(l)} ), J\{J(l)} )
  ID_to_sub_combi_infos::Vector{Tuple{Int, Int, Combination}}
end

  # _ID(combi) -> combination_index(combi)
  #ID_to_combi_index::Vector{Int}


## Display
#function Base.show(io::IO, combi::Combination{D}) where D
#  k = length(combi)
#  iszero(k) && print(io,"∅")
#  print(io,"$(join(combi.data))")
#end
#function Base.show(io::IO, ::MIME"text/plain", combi::Combination{D}) where D
#  k = length(combi)
#  if iszero(k)
#    print(io,"Combination{$D}(∅)")
#  else
#    print(io,"Combination{$D}($(join(combi.data,",")))")
#  end
#end

"""
    _ID(combi::Combination)

Internal iddentifyer of a combination of length `k`, that it is independant of
the struct parameter `D`.

This is slightly cheaper than a [`Base.hash`](@ref), and gives indices in
1:binomial(max(`combi`),k) usable to store combinations in arrays.
"""
_ID(data::Combination) = mapreduce(t -> 1<<(t-1), +, data, init=0) + 1

function _ID_to_combi(ID::Integer)
  @check ID > 0
  bits = ID-1
  k = count_ones(bits)
  combi = Vector{Int}(undef, k)
  @inbounds _ID_to_combi!(combi, bits, k)
end
function _ID_to_combi!(combi::Combination, ID::Integer)
  @check ID > 0
  bits = ID-1
  k = count_ones(bits)
  @boundscheck checkbounds(combi,k)
  @inbounds _ID_to_combi!(combi, bits, k)
end
@propagate_inbounds function _ID_to_combi!(combi, bits, k)
  curr, pos = 0, 1
  while curr < k
    if bits & 1 ≠ 0
      curr += 1
      @inbounds combi[curr] = pos
    end
    bits = bits >> 1
    pos += 1
  end
  combi
end

_ID_to_k(ID::Integer) = count_ones(ID-1)

"""
    Base.isless(a::Combination, b::Combination)

The ordering of `Combination` is right-digit to left-digit
lexicographic order, e.g. 12 < 13 < 23 < 14 < 24 < 34.

Unlike with the usual (left-digit to right-digit) lexicographic order, the
index of sorted `Combination`s of same length `k` is independant of the
dimension `D`, see [`sorted_combinations`](@ref).
"""

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
sorted_combinations(D::Int,k::Int) = sorted_combinations(Val(D),Val(k))

@generated function sorted_combinations(::Val{D},::Val{k}) where {D,k}
  iszero(k) && return :( return [ () ] )
  comp_rev_perm(tup) = ntuple(i -> D-tup[k-i+1]+1, Val(k))
  inc_perms = combinations(1:D,k) .|> (tup -> comp_rev_perm(tup)) |> reverse
  :( return $inc_perms )
end

"""
    complement(combi, Val(D))

TBW
"""
function complement(combi, ::Val{D}) where D
  k = length(combi)
  comp = Combination(undef, D-k)
  curr_perm, curr_comp = 1, 1
  for i in 1:D
    if curr_perm ≤ k && combi[curr_perm] == i
      curr_perm += 1
    else
      comp[curr_comp] = i
      curr_comp += 1
    end
  end
  comp
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
    combination_index(combi)

Linear index of `combi` amongst combinations of the same size `k`,
sorted in right-to-left lexicographic order. It depends on `k` but not on the
space dimension, see [`sorted_combinations`](@ref).
"""
@inline function combination_index(combi)
  k = length(combi)
  return sum( binomial(combi[i]-1, i) for i in 1:k; init=0) + 1
end

"""
    sub_combinations(combi::Combination{D})

Return a vector containing the `k-1` combinations `combi\\combi[i]` for 1 ≤ i ≤ `k`,
where `k=length(combi)`.
"""
function sub_combinations(combi::Combination)
  k = length(combi)
  iszero(k) && return Combination()
  sub_combis = Vector{Combination}(undef, k)
  for i in 1:k
    sub_combi = Int[combi[j + Int(j≥i)] for j in 1:(k-1)]
    sub_combis[i] = sub_combi
  end
  return sub_combis
end


#TODO
function _sub_combinations_ids!(sub_ids::AbstractVector{Tuple{Int,Int}}, combi::Combination)
  k = length(combi)
  @check length(sub_ids) >= k
  sub_combi = MVector{k-1,Int}(undef)
  for i in 1:k
    sub_combi .= ntuple(j -> combi[j + Int(j≥i)],k-1)
    sub_combi_id = combination_index(sub_combi)
    sub_ids[i] = (i, sub_combi_id)
  end
  k
end

function combination_index(combi::Combination, combi_ids::PLambdaIndices)
  id = _ID(combi)
  combi_ids.ID_to_combi_index[id]
end

function sub_combinations(combi::Combination, combi_ids::PLambdaIndices)
  id = _ID(combi)
  combi_ids.ID_to_sub_combi_infos[id]
end

function generate_ID_to_combi_index(D)
  nb_precomp_iperms = 2^D
  ID_to_combi_index = Vector{Int}(undef, nb_precomp_iperms)

  curr_indices = zero(MVector{D+1,Int})
  for ID in 1:nb_precomp_iperms
    k = _ID_to_k(ID)
    curr_indices[k+1] += 1
    ID_to_combi_index[ID] = curr_indices[k+1]
  end

  ID_to_combi_index
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
  s
end

function minor(M,I,J,::Val{k}) where k
  @check length(I) == length(J)
  @check I ⊆ axes(M)[1]
  @check J ⊆ axes(M)[2]

  T = eltype(M)
  m = MMatrix{k,k,T}(undef)
  for (i, Ii) in enumerate(take(I,k))
    for (j, Jj) in enumerate(take(J,k))
      @inbounds m[i,j] = M[Ii,Jj]
    end
  end
  det(m)
end

function all_k_minors!(m,M,Vk::Val)
  D = size(M)[1]
  Λᵏᴰ = sorted_combinations(Val(D),Vk)
  @inbounds begin
    for (i, I) in enumerate(Λᵏᴰ)
      for (j, J) in enumerate(Λᵏᴰ)
        m[i,j] = @inline minor(M,I,J,Vk)
      end
    end
  end
end

