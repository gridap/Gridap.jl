
using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Helpers
using ForwardDiff
using BenchmarkTools
using ProfileView
using StaticArrays
using Combinatorics


N = 1000
D = 2; V = Float64; order = 3

x = [rand(Point{D,Float64}) for _ in 1:N]

b = BernsteinBasisOnSimplex(Val(D),V,order)

#using Gridap.Polynomials: _upwards_de_Casteljau_indices
using Gridap.Polynomials: _downwards_de_Casteljau_indices
#using Gridap.Polynomials: _upwards_de_Casteljau_nD!# (c, λ, ::Val{K}, ::Val{D}) where {K,D}
using Gridap.Polynomials: _downwards_de_Casteljau_nD!# (c, λ, ::Val{K}, ::Val{D}) where {K,D}
using Gridap.Polynomials: _cart_to_bary # (c, λ, ::Val{K}, ::Val{D}) where {K,D}
using Gridap.Polynomials: bernstein_terms # (c, λ, ::Val{K}, ::Val{D}) where {K,D}
using Gridap.Polynomials: _sub_sub_multi_indices

function _up_de_Casteljau_nD_cost(K,D)
  count = 0
  for Ki in 1:K
    for (id,sub_ids) in _upwards_de_Casteljau_indices(Ki,D)
      for (idα⁻, d) in sub_ids
        count += 1
      end
    end
  end
  return count
end

function _dn_de_Casteljau_nD_cost(K,D)
  count = 0
  for Ki in (K-1):-1:0
    for (id,sup_ids) in _downwards_de_Casteljau_indices(Ki,D)
      for (idα⁺, d) in sup_ids
        count += 1
      end
    end
  end
  return count
end

D = 3
c = zeros(11)
K = 2
λ1 = _cart_to_bary( Point(ntuple(i -> ==(i,1) ? 1 : 0,D)), nothing)
λ2 = _cart_to_bary( Point(ntuple(i -> ==(i,2) ? 1 : 0,D)), nothing)
λ1 = ( 0., .1, 0., -1. )
λ2 = ( .0, .0, 1., -1. )
λ3 = ( 1., .0, 0., -1. )

c[1] = 1
#_upwards_de_Casteljau_nD!(c,λ3,Val(K),Val(D))

function _der_indices(K,D,λ...)
  c = zeros(binomial(K+D,D))
  c[1] = 1
  println(c)
  for (i,λi) in enumerate(λ)
    _downwards_de_Casteljau_nD!(c,λi,Val(K-i+1),Val(D),Val(K-i))
    println(c)
  end
end

#_der_indices(2,3,λ2,λ2,λ2)

function _der_indices2(c,D,λ...)
  c[1] = 1
  for (i,λi) in enumerate(λ)
    #println(i)
    _upwards_de_Casteljau_nD!(c,λi,Val(i),Val(D),Val(i))
    #println(c)
  end
  println(c)
end

#D = 2
#K = 1
#λ = [ ntuple(i -> Int(==(i,j)), D+1) for j in 1:D ]
#∇λ = [ ntuple(i -> Int(==(i,j)) - Int(==(i,D+1)), D+1) for j in 1:D ]
#c = zeros(binomial(K+D,D))
#for i in 1:D
#  print("$(λ[i])   ")
#  _der_indices2(c,D,∇λ[i])
#end

#println()
#println()
#K = 2
#c = zeros(binomial(K+D,D))
#for i in 1:D
#  for j in 1:D
#    print("$(@. λ[i]+λ[j])  ")
#    _der_indices2(c,D,∇λ[i],∇λ[j])
#  end
#end

#println()
#println()
#K = 3
#c = zeros(binomial(K+D,D))
#for i in 1:D
#  for j in 1:D
#    for k in 1:D
#      print("$(@. λ[i]+λ[j]+λ[k])  ")
#      _der_indices2(c,D,∇λ[i],∇λ[j],∇λ[k])
#    end
#  end
#end

#for term in bernstein_terms(Val(4),Val(2))
#  α = term .- 1
#  print(α)
#  for ss in _sub_sub_multi_indices(α)
#    print(" ")
#    print(ss)
#  end
#  println()
#end


module IncreasingPermutations
using Core: UnsignedMultiplicativeInverse
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Helpers
using ForwardDiff
using BenchmarkTools
using ProfileView
using StaticArrays
using Combinatorics

import Base

export IncreasingPermutation
export linear_index
export minus_one_iperms
export complement
export increasing_perms

const MAX_D = 8
const NB_IPERMS = 2^MAX_D

struct IncreasingPermutation{k,D} <: AbstractSet{Int}
  bits::UInt8

  function IncreasingPermutation{k,D}(bits::Unsigned) where {k,D}
    @check 0 ≤ D ≤ MAX_D "Only permutation of ⟦1,`D`≤$MAX_D⟧ are currently supported, got D=$D"
    @check 0 ≤ k ≤ D "0 ≤ k ≤ D is required, got k=$k and D=$D"
    mask_D = ~((~zero(bits))<<D) # mask = 0b0...01...1  with D ones
    used_bits = bits & mask_D    # keep relevant bits
    @check count_ones(used_bits)==k "permutation `bits` does not countain k=$k ones in its first D=$D bits"
    new{k,D}(bits)
  end
  function IncreasingPermutation{k,D}(perm) where {k,D}
    @check 0 ≤ D ≤ MAX_D "Only permutation of ⟦1,`D`≤$MAX_D⟧ are currently supported, got D=$D"
    @check 0 ≤ k ≤ D "0 ≤ k ≤ D is required, got k=$k and D=$D"
    @check k == length(perm) "`perm` is not of length k=$k"
    @check all(@. 1 ≤ perm ≤ D) "`perm` is not a permutation of ⟦1,`D`⟧"
    iszero(k) && return new{k,D}(0)

    bits = mapreduce(t -> 1<<(t-1), +, perm)
    @check count_ones(bits)==k "`perm` was not (strictly) increasing"
    new{k,D}(bits)
  end
end

Base.convert(::Type{IncreasingPermutation{k,D}}, x::Unsigned) where {k,D} = IncreasingPermutation{k,D}(x)

IncreasingPermutation(perm, D) = IncreasingPermutation(perm, Val(D))
IncreasingPermutation(perm, ::Val{D}) where D = IncreasingPermutation{length(perm),D}(perm)
IncreasingPermutation(perm::Unsigned, D) = IncreasingPermutation(perm, Val(D))
function IncreasingPermutation(bits::Unsigned, ::Val{D}) where D
  mask_D = ~((~zero(bits))<<D) # mask = 0b0...01...1  with D ones
  used_bits = bits & mask_D    # keep relevant bits
  IncreasingPermutation{count_ones(used_bits),D}(used_bits)
end

# Iteration Interface
Base.length(::IncreasingPermutation{k}) where k = k
Base.iterate(::IncreasingPermutation{0}) = nothing
Base.isdone(::IncreasingPermutation{k}, state) where k = state[1]==k

function Base.iterate(iperm::IncreasingPermutation{k,D}) where {k,D}
  bits = iperm.bits
  curr, pos = 0, 0
  while pos ≤ D
    pos += 1
    if bits & 1 ≠ 0
      curr += 1
      return (pos, (curr, pos))
    end
    bits = bits >> 1
  end
  @unreachable
end

function Base.iterate(iperm::IncreasingPermutation{k,D}, state) where {k,D}
  curr, pos = state
  curr == k && return nothing

  bits = iperm.bits >> pos
  while pos ≤ D
    pos += 1
    if bits & 1 ≠ 0
      curr += 1
      return (pos, (curr, pos))
    end
    bits = bits >> 1
  end
  @unreachable
end

# Indexing Interface
Base.IndexStyle(::IncreasingPermutation) = IndexLinear()
Base.keys(::IncreasingPermutation{k}) where k = Base.OneTo(k)
Base.firstindex(::IncreasingPermutation) = 1
Base.lastindex(::IncreasingPermutation{k}) where k = k

function Base.getindex(iperm::IncreasingPermutation{k}, i) where k
  @check 1 ≤ i ≤ k "Indices in IncreasingPermutations{$k} are 1:$k"

  bits = iperm.bits
  curr, pos = 0, 0
  while curr < i
    pos += 1
    if bits & 1 ≠ 0
      curr += 1
    end
    bits = bits >> 1
  end
  return Int(pos)
end

# Iterable Collection API

function Base.in(item, ipa::IncreasingPermutation)
  isint = try isinteger(item) catch _ return false end
  !isint && return false
  return ipa.bits & (1<<(Int(item)-1)) ≠ 0
end

function Base.indexin(items, iperm::IncreasingPermutation)
  indexes = similar(items, Union{Nothing, Int64})
  for i in eachindex(items)
    it = items[i]

    if it ∈ iperm
      bits_before = iperm.bits << (MAX_D-Int(it))
      indexes[i] =  count_ones(bits_before)
    else
      indexes[i] = nothing
    end
  end
  return indexes
end

Base.maximum(iperm::IncreasingPermutation{k}; init=0) where k = getindex(iperm, k)
Base.maximum(iperm::IncreasingPermutation{0}) = throw(ArgumentError())
Base.maximum(iperm::IncreasingPermutation{0}; init) = Int(init)

Base.minimum(iperm::IncreasingPermutation{k}; init=0) where k = getindex(iperm, 1)
Base.minimum(iperm::IncreasingPermutation{0}) = throw(ArgumentError())
Base.minimum(iperm::IncreasingPermutation{0}; init) = Int(init)


# Set API
function Base.union(ipa::IncreasingPermutation{ka,D}, ipb::IncreasingPermutation{kb,D}) where {ka,kb,D}
  cbits = ipa.bits | ipb.bits
  k = count_ones(cbits)
  IncreasingPermutation{k,D}(cbits)
end

function Base.intersect(ipa::IncreasingPermutation{ka,D}, ipb::IncreasingPermutation{kb,D}) where {ka,kb,D}
  cbits = ipa.bits & ipb.bits
  k = count_ones(cbits)
  IncreasingPermutation{k,D}(cbits)
end


function Base.setdiff(ipa::IncreasingPermutation{ka,D}, ipb::IncreasingPermutation{kb,D}) where {ka,kb,D}
  # TODO
end

function Base.symdiff(ipa::IncreasingPermutation{ka,D}, ipb::IncreasingPermutation{kb,D}) where {ka,kb,D}
  # TODO
end

function Base.issubset(ipa::IncreasingPermutation, ipb::IncreasingPermutation)
  return ipa.bits | ipb.bits == ipb.bits
end

function Base.issetequal(ipa::IncreasingPermutation, ipb::IncreasingPermutation)
  return ipa.bits == ipb.bits
end

function Base.isdisjoint(ipa::IncreasingPermutation, ipb::IncreasingPermutation)
  return ipa.bits & ipb.bits == zero(UInt8)
end

# Other API
"""
    complement(iperm::IncreasingPermutation{k,D}) where {k,D}

TBW
"""
function complement(iperm::IncreasingPermutation{k,D}) where {k,D}
  bits = iperm.bits
  mask_D = ~((~zero(UInt8))<<D)      # mask = 0b0...01...1  with D ones
  complement_bits = (~bits) & mask_D # reverse the D first bits
  return IncreasingPermutation{D-k,D}(complement_bits)
end

"""
    linear_index(iperm::IncreasingPermutation)

TBW
"""
function linear_index(iperm::IncreasingPermutation)
  return iperm_to_k_linearindex[iperm.bits]
end

function generate_iperm_to_k_linearindex()
  res = MVector{NB_IPERMS,UInt8}(undef)

  curr_indices = zero(MVector{MAX_D+1,UInt8})
  for bits in 0:NB_IPERMS-1
    k = count_ones(bits)
    curr_indices[k+1] += 1
    res[bits+1] = curr_indices[k+1]
  end

  return SVector(res)
end

const iperm_to_k_linearindex = generate_iperm_to_k_linearindex()


"""
    minus_one_iperms(iperm::IncreasingPermutation{k,D}) where {k,D}

TBW
"""
function minus_one_iperms(::IncreasingPermutation{0,D}) where D
  minus_one_iperms_bits = iperm_to_k_minus_one_iperms[1]
  T = IncreasingPermutation{0,D}
  minus_one_iperms = reinterpret(T, minus_one_iperms_bits)
  minus_one_iperms
end
function minus_one_iperms(iperm::IncreasingPermutation{k,D}) where {k,D}
  minus_one_iperms_bits = iperm_to_k_minus_one_iperms[iperm.bits+1]
  T = IncreasingPermutation{k-1,D}
  minus_one_iperms = reinterpret(T, minus_one_iperms_bits)
  minus_one_iperms
end

function generate_iperm_to_k_minus_one_iperms()
  res = SizedVector{NB_IPERMS,Vector{UInt8}}(undef)
  for bits in 0:NB_IPERMS-1
    res[bits+1] = generate_k_minus_one_iperms(UInt8(bits))
  end
  return SVector(res)
end

function generate_k_minus_one_iperms(bits::UInt8)
  max_bit = 1 << (MAX_D-1)
  k = count_ones(bits)
  minus_one_iperms = Vector{UInt8}(undef, k)

  bits_left = bits
  curr, pos = k, MAX_D
  while curr > 0
    if bits_left & max_bit ≠ 0
      curr -= 1
      minus_one_iperms[k-curr] = bits & ~( 1<<(pos-1) )
    end
    bits_left = bits_left << 1
    pos -= 1
  end
  minus_one_iperms
end

const iperm_to_k_minus_one_iperms = generate_iperm_to_k_minus_one_iperms()


"""
    increasing_perms(::Val{k},::Val{D})

Return a tuple (I_i) of all the increasing permutation of 1:D of length k:
1≤ I_1 < ... < I_{k} ≤ N
sorted in decreasing lexicographic order, e.g.
    {12, 13, 14, 23, 24, 34}
for k=2, D=4.
"""
@generated function increasing_perms(::Val{k},::Val{D}) where {k,D}
  iszero(k) && return :( return (0,) )
  perms = Iterators.filter(issorted, permutations(1:D,k)) .|> (perm->IncreasingPermutation{k,D}(perm))
  :( return $perms )
end
increasing_perms(k,D) = increasing_perms(Val(k),Val(D))

# Display
function Base.show(io::IO, iperm::IncreasingPermutation{k,D}) where {k,D}
  iszero(k) && print(io,"∅")
  print(io,"$(join(decode(iperm)))")
end
function Base.show(io::IO, ::MIME"text/plain", iperm::IncreasingPermutation{k,D}) where {k,D}
  if iszero(k)
    print(io,"IncreasingPermutation{$k,$D}(∅)")
  else
    print(io,"IncreasingPermutation{$k,$D}($(join(decode(iperm),",")))")
  end
end

function decode(iperm::IncreasingPermutation{k,D}) where {k,D}
  bits = iperm.bits
  curr, pos = 0, 1
  decperm = MVector{k,Int}(undef)
  while curr < k
    if bits & 1 ≠ 0
      curr += 1
      decperm[curr] = pos
    end
    bits = bits >> 1
    pos += 1
  end
  Tuple(decperm)
end

end # module IncreasingPermutations

using .IncreasingPermutations
