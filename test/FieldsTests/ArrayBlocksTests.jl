module BlocksTests

using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.Fields: ArrayBlock, MockFieldArray, MockField, BroadcastingFieldOpMap, BlockMap
using Test
using FillArrays
using LinearAlgebra
#using Gridap.ReferenceFEs

b = ArrayBlock([Int[],[1,2,3,4]],Bool[0,1])
@test length(b) == 2
@test size(b) == (2,)
@test ndims(b) == 1
@test b[1] == nothing
@test b[2] == [1,2,3,4]
b2 = ArrayBlock([ rand(2,3) for i in 1:4, j in 1:3],rand([true,false],4,3))

#@show b
#display(b)
#@show b2
#display(b2)
#display(rand(3,4))
#display(rand(3))

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
f = ArrayBlock(f_array,f_touched)
ts = fill(true,length(f_touched))
u = ArrayBlock([rand(length(f_basis)),rand(length(f_basis)),rand(length(f_basis))],ts)
uh = linear_combination(u,f)
fx = evaluate(f,x)
ft = transpose(f)
ftx = evaluate(ft,x)
uhx = evaluate(uh,x)
∇f = Broadcasting(∇)(f)
∇ft = Broadcasting(∇)(ft)
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
h = ArrayBlock(h_array,h_touched)
Jt = ∇(ϕ)
invJt = Operation(Fields.pinvJt)(Jt)
Broadcasting(Operation(⋅))(invJt,∇f)
Broadcasting(Operation(⋅))(invJt,∇ft)

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
cell_u = [ArrayBlock([rand(length(f_basis)),rand(length(f_basis)),rand(length(f_basis))],ts) for i in 1:ncells]
cell_uh = lazy_map(linear_combination,cell_u,cell_f)
cell_uhx = lazy_map(evaluate,cell_uh,cell_x)
collect(cell_uhx)
cell_∇uh = lazy_map(Broadcasting(∇),cell_uh)
cell_∇uhx = lazy_map(evaluate,cell_∇uh,cell_x)
collect(cell_∇uhx)
cell_intuh = lazy_map(integrate,cell_uh,cell_x,cell_w)
collect(cell_intuh)
cell_ϕ = fill(ϕ,ncells)
cell_ϕx = lazy_map(evaluate,cell_ϕ,cell_x)
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

# Multiple args with same size and non-zero block structure
cell_rx = lazy_map(BroadcastingFieldOpMap((i,j,k)->i-j-k),cell_fx,cell_fx,cell_fx)
collect(cell_rx)
@test cell_rx[1][1] == nothing
@test cell_rx[1][2] == -cell_fx[1][2]
@test cell_rx[1][3] == -cell_fx[1][3]

# Multiple args with same size and non-zero block structure combined with
# plain arrays
cell_rx = lazy_map(BroadcastingFieldOpMap((i,j,k)->i*j*k),cell_fx,cell_ϕx,cell_fx)
collect(cell_rx)
@test cell_rx[1][1] == nothing
@test cell_rx[1][2] == cell_fx[1][2].*cell_ϕx[1].*cell_fx[1][2]
@test cell_rx[1][3] == cell_fx[1][3].*cell_ϕx[1].*cell_fx[1][3]

nfields = 2
dv = f_basis
dq = h_basis
cell_dv_sf = Fill(dv,ncells)
cell_dq_sf = Fill(dq,ncells)
cell_dv = lazy_map(BlockMap(nfields,1),cell_dv_sf)
cell_dq = lazy_map(BlockMap(nfields,2),cell_dq_sf)
@test cell_dv[end][1] === dv
@test cell_dv[end][2] === nothing
@test cell_dq[end][1] === nothing
@test cell_dq[end][2] === dq
cell_du = lazy_map(transpose,cell_dv)
cell_dp = lazy_map(transpose,cell_dq)

cell_dv⁺ = lazy_map(BlockMap(2,1),cell_dv)
cell_dv⁻ = lazy_map(BlockMap(2,2),cell_dv)
cell_dq⁺ = lazy_map(BlockMap(2,1),cell_dq)
cell_dq⁻ = lazy_map(BlockMap(2,2),cell_dq)

cell_du⁺ = lazy_map(BlockMap((1,2),1),cell_du)
cell_du⁻ = lazy_map(BlockMap((1,2),2),cell_du)
cell_dp⁺ = lazy_map(BlockMap((1,2),1),cell_dp)
cell_dp⁻ = lazy_map(BlockMap((1,2),2),cell_dp)

cell_dvx⁺ = lazy_map(evaluate,cell_dv⁺,cell_x)
cell_dvx⁻ = lazy_map(evaluate,cell_dv⁻,cell_x)
cell_dqx⁺ = lazy_map(evaluate,cell_dq⁺,cell_x)
cell_dqx⁻ = lazy_map(evaluate,cell_dq⁻,cell_x)

cell_dux⁺ = lazy_map(evaluate,cell_du⁺,cell_x)
cell_dux⁻ = lazy_map(evaluate,cell_du⁻,cell_x)
cell_dpx⁺ = lazy_map(evaluate,cell_dp⁺,cell_x)
cell_dpx⁻ = lazy_map(evaluate,cell_dp⁻,cell_x)

cell_∇v⁺ = lazy_map(Broadcasting(∇),cell_dv⁺)
cell_∇xv⁺ = lazy_map(Broadcasting(Fields.push_∇),cell_∇v⁺,cell_ϕ)
cell_∇xvx⁺ = lazy_map(evaluate,cell_∇xv⁺,cell_x)
collect(cell_∇xvx⁺)

cell_∇u⁺ = lazy_map(Broadcasting(∇),cell_du⁺)
cell_∇xu⁺ = lazy_map(Broadcasting(Fields.push_∇),cell_∇u⁺,cell_ϕ)
cell_∇xux⁺ = lazy_map(evaluate,cell_∇xu⁺,cell_x)
collect(cell_∇xux⁺)

cell_dux_jump = lazy_map(BroadcastingFieldOpMap(-),cell_dux⁺,cell_dux⁻)
collect(cell_dux_jump)

cell_int_du⁺ = lazy_map(integrate,cell_du⁺,cell_x,cell_w)
collect(cell_int_du⁺)

cell_int_dv⁺ = lazy_map(integrate,cell_dv⁺,cell_x,cell_w)
collect(cell_int_dv⁺)

#@show cell_dux_jump[end][1,1][1,1]

cell_matp = lazy_map(BroadcastingFieldOpMap(*),cell_dqx⁺,cell_dux_jump)
cell_mat = lazy_map(IntegrationMap(),cell_matp,cell_w)

A = cell_mat[end]
b = cell_int_dv⁺[end]
c = A*b
mul!(c,A,b)

c = A*A
mul!(c,A,A)

cache = return_cache(*,A,b)
c = evaluate!(cache,*,A,b)

k = MulAddMap(2.0,1.5)
cache = return_cache(k,A,b,c)
d = evaluate!(cache,k,A,b,c)

@test d == d - d + d

@test transpose(A) != nothing

A11 = A[1,1]
A11t = transpose(A11)

@test A11*A11t != nothing
@test A*transpose(A) != nothing

I11 = ArrayBlock([[1,2,2,1],[3,1,6,2]],[true,true])
J11 = ArrayBlock([[3,2,1,5],[6,1,3,4]],[true,true])
I = ArrayBlock([I11,I11],[true,true])
J = ArrayBlock([J11,J11],[true,true])
c11 = c[1,1]

mat = zeros(ComplexF64,6,6)
vec = zeros(ComplexF64,6)
k! = AddEntriesMap(+)
k!(mat,A11,I11,J11)
k!(vec,c11,I11)
k!(mat,A,I,J)
k!(vec,c,I)
k! = TouchEntriesMap()
k!(mat,A11,I11,J11)
k!(vec,c11,I11)
k!(mat,A,I,J)
k!(vec,c,I)

#display(A)
#display(c)
#display(A[1,1])
#display(c[1])
#display(c[1][2])
end # module
