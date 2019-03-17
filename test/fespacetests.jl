using Numa, Test

##
spdims = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,spdims)
order=1
orders=order*ones(Int64,spdims)
polytope = Polytope(ones(Int64,spdims))
reffe = LagrangianRefFE(polytope,orders)
mesh = StructHexMesh(nparts)
fesp = FESpace(reffe,mesh)
@test fesp.l2giddof[end][end]==(nparts1d+1)^spdims
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
quad=Quadrature(gps)
shfscal = shfsps(reffe,quad.points)
reffe = LagrangianRefFE(polytope,orders,2)
shftens = shfsps(reffe,quad.points)
@test shftens[36,1,2,2]==shfscal[9,1,1]
##
points=quad.points
reffe = LagrangianRefFE(polytope,orders,0)
gradscal = gradshfsps(reffe,points)
reffe = LagrangianRefFE(polytope,orders,1)
gradtens = gradshfsps(reffe,points)
@test gradtens[10,1,2,2]==gradscal[1,1,2]
@test gradtens[10,3,2,2]==gradscal[1,3,2]
@test gradtens[1,3,2,1]==gradscal[1,3,2]
@test gradtens[1,3,2,2]==0.0
@test gradtens[10,3,2,1]==0.0
##
