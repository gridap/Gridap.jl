using Numa, Test

##
orders=[2,2]
polytope = Polytope([1,1])
#2*orders+1<2*gps
gps=[4,4]
quad=TensorProductQuadratureOld(gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.points)
points = quad.points
spdims = length(reffe.prebasis.polynomials)
grad  = gradient(reffe.prebasis,points)
grad = reshape(grad,size(grad,1),size(grad,2)*size(grad,3))
grad = reffe.changeofbasis*grad
grad = reshape(grad,size(grad,1),spdims,:)
##
##
gradshfs = gradshfsps(reffe,quad.points)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:prod(gps)]
elmat = sum(quad.tpweights.*elmatgp)
@test 1.0+sum(elmat)≈1.0
@test sum(elmat,dims=2).+1≈fill(0,size(elmat,1)).+1
@test elmat[1,1] ≈ 28/45
##

##
orders=[3,3,3]
polytope = Polytope([1,1,1])
#2*orders+1<2*gps
gps=[4,4,4]
quad=TensorProductQuadratureOld(gps)
reffe = LagrangianRefFE(polytope,orders)
gradshfs = gradshfsps(reffe,quad.points)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:prod(gps)]
elmat = sum(quad.tpweights.*elmatgp)
@test 1.0+sum(elmat)≈1.0
@test sum(elmat,dims=2).+1≈fill(0,size(elmat,1)).+1
##

##
polytope = Polytope([1])
orders=[1]; gps=[1]
quad=TensorProductQuadratureOld(gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.points)
elmass = [ view(shfs,:,igp)*view(shfs,:,igp)' for igp=1:prod(gps)]

gradshfs = gradshfsps(reffe,quad.points)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:prod(gps)]
ellap = sum(quad.tpweights.*elmatgp)
##

##
spdims = 2
nparts1d = 1
nparts = nparts1d*ones(Int64, spdims)
order=1
orders=order*ones(Int64, spdims)
polytope = Polytope(ones(Int64, spdims))
reffe = LagrangianRefFE(polytope, orders)
mesh = StructHexMesh(nparts)
fesp = FESpace(reffe, mesh)
u = TestFEBasis(fesp)
v = TrialFEBasis(fesp)
##
#TrialSpace(u) should allow not to compute things twice
gps1d=2
gps=gps1d*ones(Int64, spdims)
quad=TensorProductQuadratureOld(gps)
feform = v*u
sysmat = assemblesystem(feform, quad)
feform = v*u
sysmat2 = assemblesystem(feform, quad)
@test sysmat == sysmat2
@test sysmat ≈ [4/9 2/9 2/9 1/9; 2/9 4/9 1/9 2/9; 2/9 1/9 4/9 2/9; 1/9 2/9 2/9 4/9]
##
feform = u*v+∇(u)*∇(v)
sysmat = assemblesystem(feform, quad)
@test sysmat[1,1] ≈ 2/3+4/9
feop(x,y) = ∇(x)*∇(y)
sysmat = assemblesystem(feop(u,v), quad)
@test 1.0+sum(sysmat)≈1.0
@test sum(sysmat,dims=2).+1≈fill(0,size(sysmat,1)).+1
@test sysmat[1,1] ≈ 2/3
##


# Now we play with binary tree expression over BilinearForm, binary tree expressions over TrialSpaceField/TestSpaceField. Later consider functions and constants (as a type of field of 1 dim in terms of nshf). TestSpaceField could be just the "adjoint" of TrialSpaceField and thus, not repeat all operations (?). Similar as adjoint (dual here) for matrices in Julia Base.



##
# We want
#integrate(∇(u)*∇(v) + u*v, cellvalues)
#where cellvalues provide the "cells" to iterate and for each cell the integration points to be used
# I like v = TestFEBasis, u = TrialFEBasis at user level
# Next, thing how to write inside...
#
##
