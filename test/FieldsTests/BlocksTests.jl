module BlocksTests

using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: GBlock, MockFieldArray, MockField, BroadcastingFieldOpMap
using Test
using FillArrays
#using Gridap.ReferenceFEs

b = GBlock([Int[],[1,2,3,4]],Bool[0,1])
@test length(b) == 2
@test size(b) == (2,)
@test ndims(b) == 1
@test b[1] == nothing
@test b[2] == [1,2,3,4]

@test typeof(testvalue(typeof(b))) == typeof(b)

x = [Point(0.5,0.5),Point(1.0,0.0)]
w = [1.0,1.0]
#reffe = ReferenceFE(QUAD,lagrangian,Float64,1)
#f_basis = get_prebasis(reffe)
f_basis = MockFieldArray(rand(4))
f_array = Vector{typeof(f_basis)}(undef,3)
f_array[2] = f_basis
f_array[3] = f_basis
f_touched = [false,true,true]
f = GBlock(f_array,f_touched)
ts = fill(true,length(f_touched))
u = GBlock([rand(length(f_basis)),rand(length(f_basis)),rand(length(f_basis))],ts)
uh = linear_combination(u,f)
fx = evaluate(f,x)
ft = transpose(f)
ftx = evaluate(ft,x)
uhx = evaluate(uh,x)
∇f = Broadcasting(∇)(f)
∇fx = evaluate(∇f,x)
∇uh = Broadcasting(∇)(uh)
∇uhx = evaluate(∇uh,x)
intf = integrate(f,x,w)
ϕ = MockField(VectorValue(2.1,3.4))
g = Broadcasting(∘)(f,ϕ)
#gx = evaluate(g,x)
h_basis = MockFieldArray([1.0+im,3.0*im,-1.0*im,1.0*im])
h_array = Vector{typeof(h_basis)}(undef,3)
h_array[1] = h_basis
h_array[3] = h_basis
h_touched = [true,false,true]
h = GBlock(h_array,h_touched)

ncells = 10
cell_x = Fill(x,ncells)
cell_w = Fill(w,ncells)
cell_f = Fill(f,ncells)
cell_fx = lazy_map(evaluate,cell_f,cell_x)
collect(cell_fx)
cell_∇f = lazy_map(Broadcasting(∇),cell_f)
cell_∇fx = lazy_map(evaluate,cell_∇f,cell_x)
collect(cell_∇fx)
cell_ft = lazy_map(transpose,cell_f)
cell_ftx = lazy_map(evaluate,cell_ft,cell_x)
collect(cell_ftx)
cell_u = [GBlock([rand(length(f_basis)),rand(length(f_basis)),rand(length(f_basis))],ts) for i in 1:ncells]
cell_uh = lazy_map(linear_combination,cell_u,cell_f)
cell_uhx = lazy_map(evaluate,cell_uh,cell_x)
collect(cell_uhx)
cell_∇uh = lazy_map(Broadcasting(∇),cell_uh)
cell_∇uhx = lazy_map(evaluate,cell_∇uh,cell_x)
collect(cell_∇uhx)
cell_intuh = lazy_map(integrate,cell_uh,cell_x,cell_w)
collect(cell_intuh)
cell_ϕ = fill(ϕ,ncells)
cell_g = lazy_map(Broadcasting(∘),cell_f,cell_ϕ)
cell_gx = lazy_map(evaluate,cell_g,cell_x)
collect(cell_gx)

cell_h = fill(h,ncells)
cell_hx = lazy_map(evaluate,cell_h,cell_x)
cell_rx = lazy_map(BroadcastingFieldOpMap(-),cell_fx,cell_hx)
collect(cell_rx)
@test cell_rx[1][1] == -cell_hx[1][1]
@test cell_rx[1][2] == cell_fx[1][2]
@test cell_rx[1][3] == cell_fx[1][3] - cell_hx[1][3]

cell_mx =  lazy_map(BroadcastingFieldOpMap(*),cell_rx,cell_ftx)
collect(cell_mx)
@test cell_mx[end][1,1] == nothing
@test cell_mx[end][2,1] == nothing
@test cell_mx[end][3,1] == nothing

cell_mx =  lazy_map(BroadcastingFieldOpMap(*),cell_ftx,cell_rx)
collect(cell_mx)
@test cell_mx[end][1,1] == nothing
@test cell_mx[end][2,1] == nothing
@test cell_mx[end][3,1] == nothing

end # module
