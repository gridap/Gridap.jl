module MonomialBasisTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.Arrays: testvalue, Broadcasting

using Gridap.Polynomials: _q_filter, _qh_filter, _p_filter, _ph_filter

@test isHierarchical(Monomial) == true

xi = Point(2,3)
np = 5
x = fill(xi,np)

# order 0 degenerated case

order = 0
V = Float64
G = gradient_type(V,xi)
H = gradient_type(G,xi)
b = MonomialBasis(Val(2),V,order)
@test testvalue(typeof(b)) isa typeof(b)
@test get_order(b) == 0
@test get_orders(b) == (0,0)

v = V[1.0,]
g = G[(0.0, 0.0),]
h = H[(0.0, 0.0, 0.0, 0.0),]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
Hbx = repeat(permutedims(h),np)
test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

# Real-valued Q space with isotropic order

order = 1
V = Float64
G = gradient_type(V,xi)
H = gradient_type(G,xi)
b = MonomialBasis(Val(2),V,order)

v = V[1.0, 2.0, 3.0, 6.0]
g = G[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (3.0, 2.0)]
h = H[(0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0), (0.0, 1.0, 1.0, 0.0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
Hbx = repeat(permutedims(h),np)
test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

# Real-valued Q space with an isotropic order

orders = (1,3)
V = Float64
G = gradient_type(V,xi)
H = gradient_type(G,xi)
b = MonomialBasis(Val(2),V,orders)

v = V[1.0, 2.0, 3.0, 6.0, 9.0, 18.0, 27.0, 54.0]
g = G[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (3.0, 2.0), (0.0, 6.0), (9.0, 12.0), (0., 27.0), (27.0, 54.0)]
h = H[(0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0), (0.0, 1.0, 1.0, 0.0), (0.0, 0.0, 0.0, 2.0), (0.0, 6.0, 6.0, 4.0), (0.0, 0.0, 0.0, 18.0), (0.0, 27.0, 27.0, 36.0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
Hbx = repeat(permutedims(h),np)
test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

# Vector-valued Q space with isotropic order

order = 1
V = VectorValue{3,Float64}
G = gradient_type(V,xi)
H = gradient_type(G,xi)
b = MonomialBasis(Val(2),V,order)

v = V[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
      [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0],
      [3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0],
      [6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 6.0]]

g = G[[0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0],
      [0.0 0.0 0.0; 0.0 0.0 0.0], [1.0 0.0 0.0; 0.0 0.0 0.0],
      [0.0 1.0 0.0; 0.0 0.0 0.0], [0.0 0.0 1.0; 0.0 0.0 0.0],
      [0.0 0.0 0.0; 1.0 0.0 0.0], [0.0 0.0 0.0; 0.0 1.0 0.0],
      [0.0 0.0 0.0; 0.0 0.0 1.0], [3.0 0.0 0.0; 2.0 0.0 0.0],
      [0.0 3.0 0.0; 0.0 2.0 0.0], [0.0 0.0 3.0; 0.0 0.0 2.0]]

h = H[
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0),
      (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
Hbx = repeat(permutedims(h),np)
test_field_array(b,x,bx,grad=∇bx,gradgrad=Hbx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:],gradgrad=Hbx[1,:])

# Vector-valued Q space with an-isotropic order

orders = (1,2)
V = VectorValue{2,Float64}
G = gradient_type(V,xi)
b = MonomialBasis(Val(2),V,orders)

v = V[
 (1.0, 0.0), (0.0, 1.0), (2.0, 0.0), (0.0, 2.0),
 (3.0, 0.0), (0.0, 3.0), (6.0, 0.0), (0.0, 6.0),
 (9.0, 0.0), (0.0, 9.0), (18.0, 0.0), (0.0, 18.0)]

g = G[
  (0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0),
  (1.0, 0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0),
  (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 0.0, 1.0),
  (3.0, 2.0, 0.0, 0.0), (0.0, 0.0, 3.0, 2.0),
  (0.0, 6.0, 0.0, 0.0), (0.0, 0.0, 0.0, 6.0),
  (9.0, 12.0, 0.0, 0.0),(0.0, 0.0, 9.0, 12.0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])

# Real-valued P space

order = 1
V = Float64
G = gradient_type(V,xi)
filter = (e,o) -> sum(e) <= o
b = MonomialBasis(Val(2),V,order,filter)

v = V[1.0, 2.0, 3.0]
g = G[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])

# Vector-valued P space

order = 1
V = VectorValue{3,Float64}
G = gradient_type(V,xi)
filter = (e,o) -> sum(e) <= o
b = MonomialBasis(Val(2),V,order,filter)

v = V[[1.0; 0.0; 0.0], [0.0; 1.0; 0.0], [0.0; 0.0; 1.0],
      [2.0; 0.0; 0.0], [0.0; 2.0; 0.0], [0.0; 0.0; 2.0],
      [3.0; 0.0; 0.0], [0.0; 3.0; 0.0], [0.0; 0.0; 3.0]]

g = G[[0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0],
      [0.0 0.0 0.0; 0.0 0.0 0.0], [1.0 0.0 0.0; 0.0 0.0 0.0],
      [0.0 1.0 0.0; 0.0 0.0 0.0], [0.0 0.0 1.0; 0.0 0.0 0.0],
      [0.0 0.0 0.0; 1.0 0.0 0.0], [0.0 0.0 0.0; 0.0 1.0 0.0],
      [0.0 0.0 0.0; 0.0 0.0 1.0]]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])

# SymTensor-valued P space

order = 1
V = SymTensorValue{2,Float64}
G = gradient_type(V,xi)
filter = (e,o) -> sum(e) <= o
b = MonomialBasis(Val(2),V,order,filter)

v = V[(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
      (2.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.0),
      (3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0)]

g = G[(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
      (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
      (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
      (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
      (0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0),
      (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
      (0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
      (0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0),
      (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])

# SymTracelessTensor-valued P space

order = 1
V = SymTracelessTensorValue{2,Float64}
G = gradient_type(V,xi)
filter = (e,o) -> sum(e) <= o
b = MonomialBasis(Val(2),V,order,filter)

v = V[(1.0, 0.0), (0.0, 1.0),
      (2.0, 0.0), (0.0, 2.0),
      (3.0, 0.0), (0.0, 3.0)]

g = G[(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
      (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
      (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1., 0.0),
      (0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0),
      (0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.),
      (0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0)]

bx = repeat(permutedims(v),np)
∇bx = repeat(permutedims(g),np)
test_field_array(b,x,bx,grad=∇bx)
test_field_array(b,x[1],bx[1,:],grad=∇bx[1,:])


order = 1
b = MonomialBasis(Val(1),Float64,order)
@test evaluate(b,Point{1,Float64}[(0,),(1,)]) == [1.0 0.0; 1.0 1.0]

b = MonomialBasis(Val(0),VectorValue{2,Float64},order)
@test evaluate(b,Point{0,Float64}[(),()]) == VectorValue{2,Float64}[(1.0, 0.0) (0.0, 1.0); (1.0, 0.0) (0.0, 1.0)]

b = MonomialBasis(Val(0),TensorValue{2,2,Float64},order)
@test evaluate(b,Point{0,Float64}[(),()]) == TensorValue{2,2,Float64}[
  (1.0, 0.0, 0.0, 0.0) (0.0, 1.0, 0.0, 0.0) (0.0, 0.0, 1.0, 0.0) (0.0, 0.0, 0.0, 1.0);
  (1.0, 0.0, 0.0, 0.0) (0.0, 1.0, 0.0, 0.0) (0.0, 0.0, 1.0, 0.0) (0.0, 0.0, 0.0, 1.0)
]

b = MonomialBasis(Val(0),SymTensorValue{2,Float64},order)
@test evaluate(b,Point{0,Float64}[(),()]) == SymTensorValue{2,Float64}[
  (1.0, 0.0, 0.0) (0.0, 1.0, 0.0) (0.0, 0.0, 1.0);
  (1.0, 0.0, 0.0) (0.0, 1.0, 0.0) (0.0, 0.0, 1.0)
]

b = MonomialBasis(Val(0),SymTracelessTensorValue{2,Float64},order)
@test evaluate(b,Point{0,Float64}[(),()]) == SymTracelessTensorValue{2,Float64}[
  (1.0, 0.0) (0.0, 1.0); (1.0, 0.0) (0.0, 1.0)
]

order = 2

@test _q_filter( (1,2) ,order) == true
@test _q_filter( (2,0) ,order) == true
@test _q_filter( (2,2) ,order) == true
@test _q_filter( (1,1) ,order) == true
@test _q_filter( (3,1) ,order) == false

@test _qh_filter( (1,2) ,order) == true
@test _qh_filter( (2,0) ,order) == true
@test _qh_filter( (2,2) ,order) == true
@test _qh_filter( (1,1) ,order) == false
@test _qh_filter( (3,1) ,order) == false

@test _p_filter( (1,2) ,order) == false
@test _p_filter( (2,0) ,order) == true
@test _p_filter( (2,2) ,order) == false
@test _p_filter( (1,1) ,order) == true
@test _p_filter( (3,1) ,order) == false
@test _p_filter( (0,1) ,order) == true

@test _ph_filter( (1,2) ,order) == false
@test _ph_filter( (2,0) ,order) == true
@test _ph_filter( (2,2) ,order) == false
@test _ph_filter( (1,1) ,order) == true
@test _ph_filter( (3,1) ,order) == false

# value_type must be concrete even when user passes a non-concrete tensor type
for V in (TensorValue{2,2,Float64}, SymTensorValue{2,Float64},
          SkewSymTensorValue{2,Float64}, SymFourthOrderTensorValue{2,Float64})
  @test isconcretetype(value_type(MonomialBasis{2}(V, (1,1))))
end


# Higher-order gradients: hand-written kernels for N ≤ 2, forward-mode AD above. The AD
# evaluation reproduces the hand-written kernels exactly where both exist.

@testset "Higher-order gradients" begin
  for (D, k) in ((2, 3), (3, 2))
    b = MonomialBasis(Val(D), Float64, k)
    xs = [Point(ntuple(d -> 0.1d + 0.05i, D)) for i in 1:4]
    for N in 1:2
      fg = FieldGradientArray{N}(b)
      # Generic AD path (FieldArrays.jl), bypassing the hand-written kernels
      G = FieldGradientArray{N,<:AbstractArray{<:Field}}
      c = invoke(return_cache, Tuple{G,AbstractVector{<:Point}}, fg, xs)
      @test invoke(evaluate!, Tuple{Any,G,AbstractVector{<:Point}}, c, fg, xs) == evaluate(fg, xs)
    end
    for N in 3:4
      fg = FieldGradientArray{N}(b)
      r = evaluate(fg, xs)
      @test size(r) == (length(xs), length(b))
      @test eltype(r) == Fields.nth_gradient_type(Float64, first(xs), Val(N))
      # The N-th gradient is the gradient of the (N-1)-th
      rp = evaluate(Broadcasting(∇)(FieldGradientArray{N-1}(b)), xs)
      @test r ≈ rp
      # Single point
      @test evaluate(fg, first(xs)) ≈ r[1,:]
    end
  end
end

end # module
