##
using Numa, Test

using Numa.Quadratures
using Numa.Helpers
using Numa.Polytopes
using Numa.Polynomials
using Numa.Polytopes: PointInt
using Numa.FieldValues
using Numa.RefFEs
using Numa.Polynomials: TensorProductMonomialBasis

using Base.Cartesian
##

# Create dofbasis using node array for Lagrangian FEs

# Create MultivariatePolynomialBasisWithChangeOfBasis
# i.e., CanonicalBasis given DOFs

# nfacetoowndofs


##
# 1D reffe
D = 1
orders=[2]
extrusion = PointInt{D}(1)
polytope = Polytope(extrusion)
nodes = NodesArray(polytope,orders)
using Numa.RefFEs: LagrangianDOFBasis
dofsb = LagrangianDOFBasis{D,ScalarValue}(nodes.coordinates)
prebasis = TensorProductMonomialBasis{D,ScalarValue}(orders)
using Numa.Polynomials: evaluate
vals = Polynomials.evaluate(prebasis,dofsb.nodes)
@test vals == [1.0 1.0 1.0; -1.0 0.0 1.0; 1.0 0.0 1.0]
##

##
# 1D reffe
D = 2
orders=[1,1]
extrusion = PointInt{D}(1,1)
polytope = Polytope(extrusion)
nodes = NodesArray(polytope, orders)
using Numa.RefFEs: LagrangianDOFBasis
dofsb = LagrangianDOFBasis{D,VectorValue{D}}(nodes.coordinates)
prebasis = TensorProductMonomialBasis{D,VectorValue{D}}(orders)
using Numa.RefFEs: evaluate
res = RefFEs.nodeevaluate(dofsb,prebasis)
##
vals = evaluate(prebasis,dofsb.nodes)
b = Array{Float64,2}(undef,length(prebasis),length(prebasis))
nnodes = length(dofsb.nodes)
@assert nnodes*length(VectorValue{D}) == length(prebasis)
for j in 1:length(VectorValue{D})
	off = (0,nnodes*(j-1))
	for i in CartesianIndices(vals)
		ij = Tuple(i).+off
		b[ij...] = vals[i][j]
		# b[i.+off] = vals[i][j]
	end
end
b
vals[4,1][2]
function nodeevaluate(this::LagrangianDOFBasis,
	prebasis::MultivariatePolynomialBasis{D,T}) where {D,T}
	vals = evaluate(prebasis,this.nodes)
	b = Array{Float64,2}(undef,size(vals)...)
	for d in 1:D
	end
	return b
end

a = (0,4)
b = (3,1)
a.+b




function evaluate!(this::TensorProductMonomialBasis{D,T},
  points::AbstractVector{Point{D}}, v::AbstractArray{T,2}) where {D,T}
  tpcoor = i -> [ Point{1}(p[i]) for p in points]
  cooruv = [tpcoor(i) for i in 1:D]
  univals = [evaluate(this.univariatebases[i],cooruv[i]) for i in 1:D]
  cid = ntuple(i -> 1:length(this.univariatebases[i]), D)
  lent = length(T)
  cid = (cid..., 1:lent)
  cid = CartesianIndices(cid)
  for (i,j) in enumerate(cid)
    d = j[D+1]
    for k in 1:length(points)
      val = prod([ univals[i][j[i],k] for i in 1:D ])
      v[i,k] = T(ntuple(i->(i==d) ? val : 0.0, lent)...)
    end
  end
end










reffe = LagrangianRefFE(polytope,orders)
@test reffe.changeofbasis==[0.0  -0.5   0.5; 1.0   0.0  -1.0; 0.0   0.5   0.5]
##


##
# 1D reffe
D = 1
orders=[2]
extrusion = PointInt{D}(1)
polytope = Polytope(extrusion)
nodes=NodesArray(polytope,orders)
prebasis=TensorProductPolynomialBasis(orders)
reffe = LagrangianRefFE(polytope,orders)
@test reffe.changeofbasis==[0.0  -0.5   0.5; 1.0   0.0  -1.0; 0.0   0.5   0.5]
##

##
# Integration
D = 2
orders=[1,1]
extrusion = PointInt{D}(1,1)
polytope = Polytope(extrusion)
gps=(2,2)
quad=TensorProductQuadrature(orders=gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.coords)
gradshfs = gradshfsps(reffe,quad.coords)
@test sum(gradshfs) ≈ 0.0
@test sum(shfs) ≈ prod(gps)
@test size(gradshfs) == (4,4,2)
@test size(shfs) == (4,4)
##
D = 1
orders=[1]
extrusion = PointInt{D}(1)
polytope = Polytope(extrusion)
#orders/2=gps
gps=(2,)
quad=TensorProductQuadrature(orders=gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.coords)
elmatgp=[ shfs[:,igp]*shfs[:,igp]' for igp=1:prod(gps)]
elmat = sum(quad.weights.*elmatgp)
@test elmat≈[2/3 1/3; 1/3 2/3]
##

##
D=2
orders=[1,1]
extrusion = PointInt{D}(1,1)
polytope = Polytope(extrusion)
#orders/2=gps
gps=(2,2)
quad=TensorProductQuadrature(orders=gps)
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.coords)
elmatgp=[ shfs[:,igp]*shfs[:,igp]' for igp=1:prod(gps)]
elmat = sum(quad.weights.*elmatgp)
@test elmat≈[4/9 2/9 2/9 1/9; 2/9 4/9 1/9 2/9; 2/9 1/9 4/9 2/9; 1/9 2/9 2/9 4/9]
##

##
D=2
orders=[2,2]
extrusion = PointInt{D}(1,1)
polytope = Polytope(extrusion)
#orders/2=gps
# gps=(3,3)
quad=TensorProductQuadrature(orders=tuple(2*orders...))
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.coords)
numgps = length(quad.weights)
elmatgp=[ shfs[:,igp]*shfs[:,igp]' for igp=1:numgps]
elmat = sum(quad.weights.*elmatgp)
@test sum(elmat)≈4
@test elmat[1,1] ≈ 16/225
@test elmat[1,2] ≈ 8/225
##

##
D=3
orders=[2,2,2]
extrusion = PointInt{D}(1,1,1)
polytope = Polytope(extrusion)
#orders/2=gps
quad=TensorProductQuadrature(orders=tuple(2*orders...))
reffe = LagrangianRefFE(polytope,orders)
shfs = shfsps(reffe,quad.coords)
numgps = length(quad.weights)
elmatgp=[ shfs[:,igp]*shfs[:,igp]' for igp=1:prod(numgps)]
elmat = sum(quad.weights.*elmatgp)
@test sum(elmat)≈8
@test elmat[1,1] ≈ 64/3375
@test elmat[1,2] ≈ 32/3375
##

##
D=1
orders=[1]
extrusion = PointInt{D}(1)
polytope = Polytope(extrusion)
#orders/2=gps
quad=TensorProductQuadrature(orders=tuple(2*orders...))
reffe = LagrangianRefFE(polytope,orders)
gradshfs = gradshfsps(reffe,quad.coords)
numgps = length(quad.weights)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:numgps]
elmat = sum(quad.weights.*elmatgp)
@test elmat ≈ [1/2 -1/2; -1/2 1/2]
##

##
D=2
orders=[1,1]
extrusion = PointInt{D}(1,1)
polytope = Polytope(extrusion)
#2*orders+1<2*gps
quad=TensorProductQuadrature(orders=tuple(2*orders...))
reffe = LagrangianRefFE(polytope,orders)
gradshfs = gradshfsps(reffe,quad.coords)
numgps = length(quad.weights)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:numgps]
elmat = sum(quad.weights.*elmatgp)
@test 1.0+sum(elmat)≈1.0
@test sum(elmat,dims=2).+1≈fill(0,size(elmat,1)).+1
@test elmat[1,1] ≈ 2/3
##

##
D=2
orders=[2,2]
extrusion = PointInt{D}(1,1)
polytope = Polytope(extrusion)
#2*orders+1<2*gps
quad=TensorProductQuadrature(orders=tuple(2*orders...))
reffe = LagrangianRefFE(polytope,orders)
gradshfs = gradshfsps(reffe,quad.coords)
numgps = length(quad.weights)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:numgps]
elmat = sum(quad.weights.*elmatgp)
@test 1.0+sum(elmat)≈1.0
@test sum(elmat,dims=2).+1≈fill(0,size(elmat,1)).+1
@test elmat[1,1] ≈ 28/45
##

##
D=3
orders=[3,3,3]
extrusion = PointInt{D}(1,1,1)
polytope = Polytope(extrusion)
quad=TensorProductQuadrature(orders=tuple(2*orders...))
reffe = LagrangianRefFE(polytope,orders)
gradshfs = gradshfsps(reffe,quad.coords)
numgps = length(quad.weights)
elmatgp = [ view(gradshfs,:,igp,:)*view(gradshfs,:,igp,:)' for igp=1:numgps]
elmat = sum(quad.weights.*elmatgp)
@test 1.0+sum(elmat)≈1.0
@test sum(elmat,dims=2).+1≈fill(0,size(elmat,1)).+1
##

##
# Vector FE space
D = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
order=2
orders=order*ones(Int64,D)
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytope(extrusion)
reffe = LagrangianRefFE(polytope,orders)
gps=[2,2]
quad=TensorProductQuadrature(orders = tuple(2*orders...))
shfscal = shfsps(reffe,quad.coords)
reffe = LagrangianRefFE(polytope,orders,2)
shftens = shfsps(reffe,quad.coords)
@test shftens[36,1,2,2]==shfscal[9,1,1]
##
##
