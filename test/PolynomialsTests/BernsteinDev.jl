
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
using Base.Iterators: take

N = 1000
D = 2; V = Float64; order = 3

x = [rand(Point{D,Float64}) for _ in 1:N]

b = BernsteinBasisOnSimplex(Val(D),V,order)

#using Gridap.Polynomials: _upwards_de_Casteljau_indices
using Gridap.Polynomials: _downwards_de_Casteljau_indices
#using Gridap.Polynomials: _upwards_de_Casteljau_nD!# (c, λ, ::Val{K}, ::Val{D}) where {K,D}
using Gridap.Polynomials: _downwards_de_Casteljau_nD!# (c, λ, ::Val{K}, ::Val{D}) where {K,D}
using Gridap.Polynomials: _cart_to_bary # (c, λ, ::Val{K}, ::Val{D}) where {K,D}
import Gridap.Polynomials: bernstein_terms # (c, λ, ::Val{K}, ::Val{D}) where {K,D}
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

include("Combinations.jl")

using .Combinations

@inline function minor(M,I,J)
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
        m[i,j] = minor(M,I,J)
      end
    end
  end
end

function _bench(D=8,k=4)
  M = rand(SMatrix{D,D,Float64})
  nb_perms = binomial(D,k)
  m = MMatrix{nb_perms,nb_perms,Float64}(undef)
  Vk = Val(k)
  @btime all_k_minors!($m,$M,$(Vk))
end

bernstein_terms(k,D) = bernstein_terms(Val(k),Val(D))

function supp(α)
  s = Int[]
  for (i,αi) in enumerate(α)
    if αi > 0
      push!(s, i)
    end
  end
  Tuple(s)
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

r,k,D = 2, 2, 3

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

function P_bubles(;r=2,k=2,D=3)
  for (d, F, dF_bubbles) in PΛ_bubble_indices(r,k,D)
    s = Set()
    println("d = $d, F=$F, F*=$(complement(F))")
    for (i, α, J) in dF_bubbles
      print("i=$i, α=$α, J=$J")
      println("  REDUNDANT")
      #if α in s
      #  println("  REDUNDANT")
      #else
      #  println()
      #  push!(s,α)
      #end
    end
    println()
  end
end

function Pm_bubles(;r=2,k=2,D=3)
  for (d, F, dF_bubbles) in P⁻Λ_bubble_indices(r,k,D)
    println("d = $d, F=$F, F*=$(complement(F))")
    for (i, α, J) in dF_bubbles
      println("i=$i, α=$α, J=$J")
      #for (l,J_l) in enumerate(sub_combinations(J))
        #println("sgn=$(-(-1)^l), J[l]=$(J[l]), J\\l=$(J_l)")
      #end
    end
    println()
  end
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
