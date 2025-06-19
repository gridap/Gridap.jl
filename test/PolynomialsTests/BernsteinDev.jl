
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

#using Gridap.Polynomials: _upwards_de_Casteljau_indices
using Gridap.Polynomials: _downwards_de_Casteljau_indices
#using Gridap.Polynomials: _upwards_de_Casteljau_nD!# (c, λ, ::Val{K}, ::Val{D}) where {K,D}
using Gridap.Polynomials: _downwards_de_Casteljau_nD!# (c, λ, ::Val{K}, ::Val{D}) where {K,D}
using Gridap.Polynomials: _cart_to_bary # (c, λ, ::Val{K}, ::Val{D}) where {K,D}
import Gridap.Polynomials: bernstein_terms # (c, λ, ::Val{K}, ::Val{D}) where {K,D}
using Gridap.Polynomials: _sub_sub_multi_indices

include("Combinations.jl")

using .Combinations

bernstein_terms(k,D) = bernstein_terms(Val(k),Val(D))

function all_k_minors!(m,M,Vk::Val)
  D = size(M)[1]
  Λᵏᴰ = _sorted_combinations(Val(D),Vk)
  @inbounds begin
    for (i, I) in enumerate(Λᵏᴰ)
      for (j, J) in enumerate(Λᵏᴰ)
        m[i,j] = @inline minor(M,I,J,Vk)
      end
    end
  end
end

"""
    _ID(combi::Combination)

Internal iddentifyer of a combination of length `k`, that it is independant of
the struct parameter `D`.

This is slightly cheaper than a [`Base.hash`](@ref), and gives indices in
1:binomial(max(`combi`),k) usable to store combinations in arrays.
"""
_ID(combi::Combination) = mapreduce(t -> 1<<(t-1), +, combi, init=0) + 1

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

function BernsteinTerm_to_Combination(α)
  k = length(α)
  iszero(k) && return []
  r = sum(α)
  v = zeros(Int,r)
  i_v = 1
  i = 1
  for α_i in α
    for _ in 1:α_i
      v[i_v] = i
      i += 1
      i_v += 1
    end
    i += 1
  end
  #println(r,k,r+k-1,v)
  return v
end

function sub_combinations(combi)
  k = length(combi)
  iszero(k) && return Vector{Int}[]
  sub_combis = Vector{Vector{Int}}(undef, k)
  for i in 1:k
    sub_combi = Int[combi[j + Int(j≥i)] for j in 1:(k-1)]
    sub_combis[i] = sub_combi
  end
  return sub_combis
end
