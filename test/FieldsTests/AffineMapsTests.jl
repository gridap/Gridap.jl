module AffineMapsTests

using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using FillArrays
using Test

origin = Point(1,1)
g = TensorValue(2,0,0,2)

h = AffineField(g,origin)
@test isa(∇(h),ConstantField)
@test isa(Broadcasting(∇)(h),ConstantField)

hinv = inverse_map(h)

x1 = Point(0,0)
x2 = Point(1,1)
x3 = Point(2,1)
x = [x1,x2,x3]
hx = h(x)
∇hx = ∇(h)(x)

r = Point{2,Int}[(1, 1), (3, 3), (5, 3)]
∇r = TensorValue{2,2,Int,4}[(2, 0, 0, 2), (2, 0, 0, 2), (2, 0, 0, 2)]
test_field(h,x,r,grad=∇r)
test_field(h,x1,x1⋅g+origin,grad=g)
test_field(hinv,r,x,grad=inv.(∇r))

ncells = 10
cell_to_h = fill(h,ncells)
cell_to_∇h = lazy_map(Broadcasting(∇),cell_to_h)
cell_to_x = Fill(x,ncells)

cell_to_hx = lazy_map(evaluate,cell_to_h,cell_to_x)
cell_to_∇hx = lazy_map(evaluate,cell_to_∇h,cell_to_x)

test_array(cell_to_hx,fill(hx,ncells))
test_array(cell_to_∇hx,fill(∇hx,ncells))

T = AffineField{3,3,Int}
@test isa(zero(T),T)

#display(cell_to_hx)
#display(cell_to_∇hx)
#print_op_tree(cell_to_hx)
#print_op_tree(cell_to_∇hx)

#using BenchmarkTools
#
#c = return_cache(h,x)
#@btime evaluate!($c,$h,$x)
#
#c = return_cache(h,x1)
#@btime evaluate!($c,$h,$x1)


# Pushforward of the N-th gradient: push_∇ and push_∇∇ at N = 1, 2

Jt_inv = TensorValue(0.5, 0.1, -0.2, 2.0)
g1 = VectorValue(1.0, -3.0)
@test PushNthGradient{1}()(g1, Jt_inv) ≈ Jt_inv ⋅ g1
g2 = TensorValue(1.0, 2.0, 2.0, -1.0)
@test PushNthGradient{2}()(g2, Jt_inv) ≈ push_∇∇(g2, Jt_inv)

end # module
