module LazyFieldArraysTests

using Test
using Gridap.Arrays
using Gridap.Fields
using Gridap.TensorValues
using FillArrays
using Gridap.Helpers

# using BenchmarkTools
using Test


np = 4
p1 = Point(1.0,2.0)
p2 = Point(2.0,1.0)
p3 = Point(4.0,3.0)
p4 = Point(6.0,1.0)

p = p1

np = 4
p1 = Point(1.0,2.0)
p2 = Point(2.0,1.0)
p3 = Point(4.0,3.0)
p4 = Point(6.0,1.0)

p = p1
np = 4
x = [p1,p2,p3,p4]

f1 = GenericField(x->1*x)
f2 = GenericField(x->2*x)
f3 = GenericField(x->3*x)
f4 = GenericField(x->4*x)
f5 = GenericField(x->5*x)

h(x) = 2.0*x
f = GenericField(h)

q(x) = Point(sqrt.(x.data))
g = GenericField(q)

nf = 2
# basis = [f1, f2, f3, f4, f5]
basis = fill(f,nf)

field = g

na = 3

x_a = fill(x,na)
basis_a = fill(basis,na)
field_a = fill(field,na)

v = rand(nf)
v_a = fill(v,na)

vm = rand(nf,nf)
vm_a = fill(vm,na)


# Evaluate field

res_field_a = lazy_map(evaluate,field_a,x_a)

c_r = array_cache(res_field_a)
# @btime getindex!($c_r,$res_field_a,1);

# Evaluate basis

res_basis_a = lazy_map(evaluate,basis_a,x_a)

c_r = array_cache(res_basis_a)
# @btime getindex!($c_r,$res_basis_a,1);

cb = return_cache(basis,x)
# @btime evaluate!(cb,basis,x)

# Transpose basis

tbasis_a = lazy_map(transpose,basis_a)
@test all([ (tbasis_a[i] == transpose(basis_a[i])) for i in 1:na])

res_tbasis_a = lazy_map(evaluate,tbasis_a,x_a)
for j in 1:na
  @test all([res_tbasis_a[j][i,1,:] == res_basis_a[j][i,:] for i in 1:np])
end

c_r = array_cache(res_tbasis_a)
# @btime getindex!($c_r,$res_tbasis_a,1);

# Broadcast operation basis basis

op = +
brbasis_a = lazy_map(Broadcasting(Operation(op)),basis_a,basis_a)
res_brbasis_a = lazy_map(evaluate,brbasis_a,x_a)

c_r = array_cache(res_brbasis_a)
# @btime getindex!($c_r,$res_brbasis_a,$1);

# Broadcast operation basis field

op = +
brbasis_a = lazy_map(Broadcasting(Operation(op)),basis_a,field_a)

res_brbasisfield_a = lazy_map(evaluate,brbasis_a,x_a)
for j in 1:na
  for i in 1:np
    @test res_brbasisfield_a[j][i,:] == res_basis_a[j][i,:] .+ res_field_a[j][i]
  end
end

c_r = array_cache(res_brbasisfield_a)
# @btime getindex!($c_r,$res_brbasisfield_a,1);

# Linear Combination basis values (vector values)

vxbasis_a = lazy_map(linear_combination,v_a,basis_a)
res_vxbasis_a = lazy_map(evaluate,vxbasis_a,x_a)

for j in 1:na
  for i in 1:np
    res_vxbasis_a[j][i] == transpose(res_basis_a[j][i,:])*v_a[j]
  end
end

c_r = array_cache(res_vxbasis_a)
# @btime getindex!(c_r,res_vxbasis_a,1)

# Linear Combination basis values (matrix values)

vmxbasis_a = lazy_map(linear_combination,vm_a,basis_a)
res_vmxbasis_a = lazy_map(evaluate,vmxbasis_a,x_a)

for j in 1:na
  for i in 1:np
    res_vmxbasis_a[j][i] == transpose(res_basis_a[j][i,:])*vm_a[j]
  end
end

c_r = array_cache(res_vmxbasis_a)
# @btime getindex!(c_r,res_vmxbasis_a,1)

# basis*transpose(basis)

basisxtbasis_a = lazy_map(Broadcasting(Operation(⋅)),basis_a,tbasis_a)
res_basisxtbasis_a = lazy_map(evaluate,basisxtbasis_a,x_a)
size(res_basisxtbasis_a[1])
for j in 1:na
  for i in 1:np
    @test res_basisxtbasis_a[j][i,:,:] == broadcast(⋅,res_basis_a[j][i,:],res_tbasis_a[j][i,:,:])
  end
end

c_r = array_cache(res_basisxtbasis_a)
# @btime getindex!($c_r,$res_basisxtbasis_a,1);

# composition basis(field)

op = ∘
compfield_a = lazy_map(op,field_a,field_a)

res_field_a = lazy_map(evaluate,field_a,x_a)


res_compfield_a = lazy_map(evaluate,compfield_a,x_a)

@test all([
  res_compfield_a[i]  ==   lazy_map(evaluate,field_a,res_field_a)[i]
  for i in 1:na])

c_r = array_cache(res_compfield_a)
# @btime getindex!($c_r,$res_compfield_a,1);

# composition basis(basis, field)

op = ∘
compbasisfield_a = lazy_map(Broadcasting(op),basis_a,field_a)

res_compbasisfield_a = lazy_map(evaluate,compbasisfield_a,x_a)
@test all([ res_compbasisfield_a[i] == lazy_map(evaluate,basis_a,res_field_a)[i] for i in 1:na])

c_r = array_cache(res_compbasisfield_a)
# @btime getindex!($c_r,$res_compbasisfield_a,1);


## transpose(basis)*basis
#
#tbasisxbasis_a = lazy_map(*,tbasis_a,basis_a)
#res_tbasisxbasis_a = lazy_map(evaluate,tbasisxbasis_a,x_a)
#for j in 1:na
#  for i in 1:np
#    res_tbasisxbasis_a[j][i] == res_tbasis_a[j][i,1,:] ⋅ res_basis_a[j][i,:]
#  end
#end
#
#c_r = array_cache(res_tbasisxbasis_a)
## @btime getindex!($c_r,$res_tbasisxbasis_a,1);

# integration (field)

w = rand(np)
w_a = fill(w,na)
jac = ConstantField(TensorValue(4.0,0.0,0.0,4.0))
jac_a = fill(jac,na)
res_jac_a = lazy_map(evaluate,jac_a,x_a)

field_a
integrate_a = lazy_map(integrate,field_a,x_a,w_a,jac_a)


# integration (basis)

# w_a = v_a
jac = ConstantField(TensorValue(4.0,0.0,0.0,4.0))
jac_a = fill(jac,na)
res_jac_a = lazy_map(evaluate,jac_a,x_a)

integrate_a = lazy_map(integrate,basis_a,x_a,w_a,jac_a)

# Gradient

∇field_a = lazy_map(gradient,field_a)
res_∇field_a = lazy_map(evaluate,∇field_a,x_a)
@test res_∇field_a[1][1] == TensorValue(0.5, 0.0, 0.0, 0.35355339059327373)

∇basis_a = lazy_map(Broadcasting(∇),basis_a)
res_∇basis_a = lazy_map(evaluate,∇basis_a,x_a)
@test all(res_∇basis_a[1] .== TensorValue(2.0, 0.0, 0.0, 2.0))

end # module
