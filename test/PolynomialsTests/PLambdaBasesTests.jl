module PLambdaBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using ForwardDiff
using StaticArrays
using BenchmarkTools

import Gridap.Arrays.evaluate!
T = Float64


#using Gridap.Polynomials: PΛ_bubble_indices, P⁻Λ_bubble_indices, complement, sorted_combinations, _hat_Ψ
#export _bench
#function _bench(D=8,k=4)
#  M = rand(SMatrix{D,D,Float64})
#  nb_perms = binomial(D,k)
#  m = MMatrix{nb_perms,nb_perms,Float64}(undef)
#  Vk = Val(k)
#  @btime all_k_minors!($m,$M,$(Vk))
#end
#
#export P_bubles
#function P_bubles(;r=2,k=2,D=3)
#  for (d, F, dF_bubbles) in PΛ_bubble_indices(r,k,D)
#    #s = Set()
#    println("d = $d, F=$F, F*=$(complement(F))")
#    for (i, α, J) in dF_bubbles
#      println("i=$i, α=$α, J=$J")
#      for I in sorted_combinations(D,k)
#        M = MMatrix{k,k,Float64}(undef)
#        for (i, Ii) in enumerate(I)
#          for (j, Jj) in enumerate(J)
#            @inbounds M[i,j] = Int((Ii+1)==Jj) - Int(Jj==1)
#          end
#        end
#        s = Int(isone(J[1]))
#        n = count(i-> (J[i]-1)∉I, (s+1):length(J))
#        p = _findfirst_or_zero(j-> (I[j]+1)∉J, 1,length(J))
#        m = _findfirst_or_zero(i-> (J[i]-1)∉I, (s+1),length(J))
#        q = _findfirst_or_zero(j-> (I[j]+s)∉J, (p+1),length(J))
#        println("\tI=$I, J=$J, M_IJ=$M, Ψ[I,J] = $(_hat_Ψ(r,α,F,I,J,Float64))")
#        println("\ts=$s, n=$n, p=$p, m=$m, q=$q")
#        println()
#      end
#      #if α in s
#      #  println("  REDUNDANT")
#      #else
#      #  println()
#      #  push!(s,α)
#      #end
#    end
#    println()
#  end
#end
#function _findfirst_or_zero(pred, start, endd)
#  r = findfirst(pred,start:endd)
#  return isnothing(r) ? 0 : r+start-1
#end
#
#export Pm_bubles
#function Pm_bubles(;r=2,k=2,D=3)
#  for (d, F, dF_bubbles) in P⁻Λ_bubble_indices(r,k,D)
#    println("d = $d, F=$F, F*=$(complement(F))")
#    for (i, α, J) in dF_bubbles
#      println("i=$i, α=$α, J=$J")
#      #for (l,J_l) in enumerate(sub_combinations(J))
#        #println("sgn=$(-(-1)^l), J[l]=$(J[l]), J\\l=$(J_l)")
#      #end
#    end
#    println()
#  end
#end

# 0D                                           0D #
D = 0
vertices = (Point{D,T}(),)
x = [vertices[1]]
x1 = x[1]

k = 0
r = 1
b = PLambdaBasis(Val(D),T,r,k)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

r = 2
b = PLambdaBasis(Val(D),T,r,k)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

r = 3
b = PLambdaBasis(Val(D),T,r,k)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

r = 4
b = PLambdaBasis(Val(D),T,r,k)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

# 1D                                           1D #
D = 1
Pt = Point{D,T}
vertices = (Pt(.5),Pt(1.))
x = [xi for xi in vertices]
x1 = x[1]

r = 4
k = 0
b = PLambdaBasis(Val(D),T,r,k)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

k = 1
b = PLambdaBasis(Val(D),T,r,k)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
evaluate(b,x)
evaluate(Broadcasting(∇)(b),x)
evaluate(Broadcasting(∇∇)(b),x)

# 2D                                           2D #
D = 2
Pt = Point{D,T}
vertices = (Pt(0., 0.5),Pt(1.,0),Pt(.5,1.))
x = [xi for xi in vertices]
x1 = x[1]

r = 4
k = 0
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

k = 1
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

k = 2
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)


# 3D                                           3D #
D = 3
Pt = Point{D,T}
vertices = (Pt(0., 0., 0.5),Pt(1.,0,0),Pt(0,.5,0),Pt(0,.5,.5))
x = [xi for xi in vertices]
x1 = x[1]

r = 4
k = 0
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

k = 1
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

k = 2
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)


k = 3
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

# 4D                                           4D #
D = 4
Pt = Point{D,T}
vertices = (Pt(0.,0.,0.,0.),Pt(0.,0.,0.,0.5),Pt(0.,1.,0,0),Pt(0.,0,.5,0),Pt(.5,1,1,1))
x = [xi for xi in vertices]
x1 = x[1]

r = 4
k = 0
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

k = 1
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

k = 2
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)


b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

k = 3

b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)


b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)

k = 4
b = PLambdaBasis(Val(D),T,r,k)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)


b = PLambdaBasis(Val(D),T,r,k,vertices)
bx  = evaluate(b,x)
Gbx = evaluate(Broadcasting(∇)(b),x)
Hbx = evaluate(Broadcasting(∇∇)(b),x)


end # module

#using .PLambdaBasisTests
#
#function α_to_I(::NTuple{0,Int})
#  return Combination{0,0}()
#end
#function α_to_I(α::NTuple{k,Int}) where k
#  r = sum(α)
#  v = zero(MVector{r,Int})
#  i_v = 1
#  i = 1
#  for α_i in α
#    for _ in 1:α_i
#      v[i_v] = i
#      i += 1
#      i_v += 1
#    end
#    i += 1
#  end
#  #println(r,k,r+k-1,v)
#  return Gridap.Polynomials.Combination{r,k+r-1}(v)
#end
