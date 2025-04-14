#######################
# PLambdaIndices type #
#######################

# just show some `PLambdaIndices` in REPL if lazy to read
# All (sub)types "XX" below should be understood as "Indices_For_XX"

# A combination iddentyfing a face or form component
# a.k.a an increasing set of indices in a range ⟦ 1, D/N/... ⟧
#     F / J / I   = 1 ≤ F1 < ... < Fd ≤ D
const Combination = Vector{Int}

# One bubble function:
#               ω_w  =      (  w,             α, α_id,           J,   sub_J_ids,  sup_α_ids )
const BubbleFunction = Tuple{Int, BernsteinTerm,  Int, Combination, Vector{Int}, Vector{Int}}

# One bubble space associated to the d-dimensional face F
# bubble_d_F =      (          F,    F_bubble_functions )
const Bubble = Tuple{Combination, Vector{BubbleFunction}}

# Indices for the components of a basis k-form or its vector calculus proxy,
#                      (I_id,           I, sgnIcomp)
const Component = Tuple{ Int, Combination,      Int}
# where (I_id, I, sgnIcomp) =
#     (combination_index(I), I, 1)                               for k-forms components,
#  or (combination_index(complement(I)), I, combination_sign(I)) for vector proxy
# for I in sorted_combinations(D,k)

struct PLambdaIndices
  # An objectid is stored here to keep track of what's inside and perform a
  # sanity check when reusing the indices elsewhere:
  #     identity = objectid( (r,k,D,:P,false) )
  # for a PᵣΛᵏ(△ᴰ) basis indices with vector-proxied components.
  identity::UInt

  # bubble indices given by P(m)Λ_bubbles (can be filtered to select bubbles)
  bubbles::Vector{Bubble}

  # Components "I" of a basis polynomial (form)
  components::Vector{Component}

  # Components "I" of the exterior derivative of a basis polynomial form
  # ext_deriv_components::Vector{Component}
end

# Display
function Base.show(io::IO, ::MIME"text/plain", indices::PLambdaIndices)
  if isempty(indices.bubbles)
    print(io,"Empty PΛ basis indices")
    return
  end

  α = indices.bubbles[1][2][1][2]
  J = indices.bubbles[1][2][1][4]
  I = indices.components[1][2]

  r = sum(α)+length(J)-length(I)
  k = length(I)
  D = length(α)-length(J)+length(I)

  println(io,"PᵣΛᵏ(△ᴰ) basis indices, r=$r k=$k D=$D")
  for (F, F_bubble) in indices.bubbles
    isempty(F_bubble) && continue

    println(io,"Bubble on face F=$(join(F))")
    println(io,"\tw \tα \tα_id \tJ \tsub_J_ids \tsup_α_ids")
    for (w, α, α_id, J, sub_J_ids, sup_α_ids) in F_bubble
      println(io,"\t$w \t$(join(α)) \t$α_id \t$(join(J)) \t$(join(sub_J_ids,",")) \t\t$(join(sup_α_ids,","))")
    end
  end

  println()
  println(io,"Basis polynomial components")
  println(io,"\tI_id\t I\t I_sgn")
  for (I_id, I, I_sgn) in indices.components
    println(io,"\t$I_id\t $(join(I))\t $I_sgn")
  end
end

"""
    _ID(combi::Combination)

Internal iddentifyer of a combination of length `k`, that it is independant of
the struct parameter `D`.

This is slightly cheaper than a [`Base.hash`](@ref), and gives indices in
1:binomial(max(`combi`),k) usable to store combinations in arrays.
"""
#_ID(combi::Combination) = mapreduce(t -> 1<<(t-1), +, combi, init=0) + 1
#
#function _ID_to_combi(ID::Integer)
#  @check ID > 0
#  bits = ID-1
#  k = count_ones(bits)
#  combi = Vector{Int}(undef, k)
#  @inbounds _ID_to_combi!(combi, bits, k)
#end
#function _ID_to_combi!(combi::Combination, ID::Integer)
#  @check ID > 0
#  bits = ID-1
#  k = count_ones(bits)
#  @boundscheck checkbounds(combi,k)
#  @inbounds _ID_to_combi!(combi, bits, k)
#end
#@propagate_inbounds function _ID_to_combi!(combi, bits, k)
#  curr, pos = 0, 1
#  while curr < k
#    if bits & 1 ≠ 0
#      curr += 1
#      @inbounds combi[curr] = pos
#    end
#    bits = bits >> 1
#    pos += 1
#  end
#  combi
#end
#
#_ID_to_k(ID::Integer) = count_ones(ID-1)

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
  :( return $(Tuple(inc_perms)) )
end

"""
    complement(combi, D)

TBW
"""
function complement(combi, D)
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

#"""
#    sub_combinations(combi::Combination{D})
#
#Return a vector containing the `k-1` combinations `combi\\combi[i]` for 1 ≤ i ≤ `k`,
#where `k=length(combi)`.
#"""
#function sub_combinations(combi::Combination)
#  k = length(combi)
#  iszero(k) && return Combination()
#  sub_combis = Vector{Combination}(undef, k)
#  for i in 1:k
#    sub_combi = Int[combi[j + Int(j≥i)] for j in 1:(k-1)]
#    sub_combis[i] = sub_combi
#  end
#  return sub_combis
#end


# TODO
function _sub_combinations_ids(combi::Combination)
  k = length(combi)
  sub_combi = MVector{k-1,Int}(undef)
  sub_combi_ids = Vector{Int}(undef, k)
  for i in 1:k
    sub_combi .= ntuple(j -> combi[j + Int(j≥i)],k-1)
    sub_combi_id = combination_index(sub_combi)
    sub_combi_ids[i] = sub_combi_id
  end
  sub_combi_ids
end

# TODO
function _sup_multi_indices(α::BernsteinTerm)
  N = length(α)
  sup_α = MVector{N,Int}(undef)
  sup_α_ids = Vector{Int}(undef, N)
  for i in 1:N
    sup_α .=  ntuple(k -> α[k]+Int(k==i), N)
    sup_α_id = bernstein_term_id(sup_α)
    sup_α_ids[i] = sup_α_id
  end
  return sup_α_ids
end

"""
    support(α)

TBW
"""
function support(α::BernsteinTerm)
  s = Int[]
  for (i,αi) in enumerate(α)
    if αi > 0
      push!(s, i)
    end
  end
  s
end

# TODO
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

