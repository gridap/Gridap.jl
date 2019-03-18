using Numa, Test

##
# 1D reffe
orders=[2]
polytope = Polytope([1])
reffe = LagrangianRefFE(polytope,orders)
@test reffe.changeofbasis==[0.0  -0.5   0.5; 1.0   0.0  -1.0; 0.0   0.5   0.5]
##

##
# Integration
orders=[1,1]
polytope = Polytope([1,1])
gps=[2,2]
quad=TensorProductQuadrature(gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.points)
gradshfs = gradshfsps(reffe,quad.points)
@test sum(gradshfs) ≈ 0.0
@test sum(shfs) ≈ prod(gps)
@test size(gradshfs) == (4,4,2)
@test size(shfs) == (4,4)
##
orders=[1]
polytope = Polytope([1])
#orders/2=gps
gps=[2]
quad=TensorProductQuadrature(gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.points)
elmatgp=[ shfs[:,igp]*shfs[:,igp]' for igp=1:prod(gps)]
elmat = sum(quad.tpweights.*elmatgp)
@test elmat≈[2/3 1/3; 1/3 2/3]
##

##
orders=[1,1]
polytope = Polytope([1,1])
#orders/2=gps
gps=[2,2]
quad=TensorProductQuadrature(gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.points)
elmatgp=[ shfs[:,igp]*shfs[:,igp]' for igp=1:prod(gps)]
elmat = sum(quad.tpweights.*elmatgp)
@test elmat≈[4/9 2/9 2/9 1/9; 2/9 4/9 1/9 2/9; 2/9 1/9 4/9 2/9; 1/9 2/9 2/9 4/9]
##

##
orders=[2,2]
polytope = Polytope([1,1])
#orders/2=gps
gps=[3,3]
quad=TensorProductQuadrature(gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.points)
elmatgp=[ shfs[:,igp]*shfs[:,igp]' for igp=1:prod(gps)]
elmat = sum(quad.tpweights.*elmatgp)
@test sum(elmat)≈4
@test elmat[1,1] ≈ 16/225
@test elmat[1,2] ≈ 8/225
##

##
orders=[2,2,2]
polytope = Polytope([1,1,1])
#orders/2=gps
gps=[6,6,6]
quad=TensorProductQuadrature(gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.points)
elmatgp=[ shfs[:,igp]*shfs[:,igp]' for igp=1:prod(gps)]
elmat = sum(quad.tpweights.*elmatgp)
@test sum(elmat)≈8
@test elmat[1,1] ≈ 64/3375
@test elmat[1,2] ≈ 32/3375
##

##
orders=[1]
polytope = Polytope([1])
#orders/2=gps
gps=[4]
quad=TensorProductQuadrature(gps)
reffe = LagrangianRefFE(polytope,orders)
gradshfs = gradshfsps(reffe,quad.points)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:prod(gps)]
elmat = sum(quad.tpweights.*elmatgp)
@test elmat ≈ [1/2 -1/2; -1/2 1/2]
##

##
orders=[1,1]
polytope = Polytope([1,1])
#2*orders+1<2*gps
gps=[2,2]
quad=TensorProductQuadrature(gps)
reffe = LagrangianRefFE(polytope,orders)
gradshfs = gradshfsps(reffe,quad.points)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:prod(gps)]
elmat = sum(quad.tpweights.*elmatgp)
@test 1.0+sum(elmat)≈1.0
@test sum(elmat,dims=2).+1≈fill(0,size(elmat,1)).+1
@test elmat[1,1] ≈ 2/3
##

##
orders=[2,2]
polytope = Polytope([1,1])
#2*orders+1<2*gps
gps=[4,4]
quad=TensorProductQuadrature(gps)
reffe = LagrangianRefFE(polytope,orders)
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
quad=TensorProductQuadrature(gps)
reffe = LagrangianRefFE(polytope,orders)
gradshfs = gradshfsps(reffe,quad.points)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:prod(gps)]
elmat = sum(quad.tpweights.*elmatgp)
@test 1.0+sum(elmat)≈1.0
@test sum(elmat,dims=2).+1≈fill(0,size(elmat,1)).+1
##

##
# Vector FE space
spdims = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,spdims)
order=2
orders=order*ones(Int64,spdims)
polytope = Polytope(ones(Int64,spdims))
reffe = LagrangianRefFE(polytope,orders)
gps=[2,2]
quad=TensorProductQuadrature(gps)
shfscal = shfsps(reffe,quad.points)
reffe = LagrangianRefFE(polytope,orders,2)
shftens = shfsps(reffe,quad.points)
@test shftens[36,1,2,2]==shfscal[9,1,1]
##
##
