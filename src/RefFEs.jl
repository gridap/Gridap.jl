using Numa.Polytopes
using Numa.Polynomials

export LagrangianRefFE
export shfsps, gradshfsps

"""
Abstract Reference Finite Element
"""
abstract type RefFE end

"""
Reference Finite Element a la Ciarlet, i.e., it relies on a local function (polynomial) space, an array of nodes (DOFs), and a polytope (cell topology). The rank of the approximating field can be arbitrary. The current implementation relies on the prebasis (e.g., monomial basis of polynomials) and a change-of-basis (using the node array) to generate the canonical basis, i.e., the shape functions.
"""
struct LagrangianRefFE{D} <: RefFE
	polytope::Polytope{D}
	prebasis::TensorProductPolynomialBasis
	nodes::NodesArray{D}
	changeofbasis::Array{Float64,2}
	#dofs::Array{Function,1} (Other way to express it... as a function?)
	# fieldrank:: scalar, vector, tensor (or better field rank)
	# conformity, continuity (here or in the DOF numbering algorithm?)
	# dofsnface, dofsclosurenface
	dofsnface::Array{Array{Int64},1}
	rank::Int64
	numdof::Int64
end

function LagrangianRefFE(polytope::Polytope{D}, orders::Array{Int64,1},
	rank::Int64=0) where D
	nodes=NodesArray(polytope,orders)
	prebasis=TensorProductPolynomialBasis(orders)
	dofspb=prebasis(nodes.coordinates)
	# dofspbs = dofbasis(prebasis)
	changeofbasis=inv(dofspb)
	dofsnface=nodes.nfacenodes
	numdof = size(changeofbasis,1)*D^rank
	LagrangianRefFE{D}(polytope, prebasis, nodes, changeofbasis, dofsnface, rank, numdof)
end

"""
Compute shape functions in a set of points
"""
function shfsps(reffe::LagrangianRefFE{D},points) where D
	scalm=reffe.changeofbasis*reffe.prebasis(points)
	nshf=size(reffe.changeofbasis,1)
	rgshfdim=[(i-1)*nshf+1:nshf*i for i=1:D^(reffe.rank)]
	aux1=D^reffe.rank
	shfss=zeros(Float64,rgshfdim[end][end],size(scalm,2),aux1)
	for i=1:length(rgshfdim)
		shfss[rgshfdim[i],:,i]=scalm
	end
	aux2=D*ones(Int64,reffe.rank)
	shfss=reshape(shfss,size(shfss,1),size(shfss,2),aux2...)
	return shfss
end

"""
Compute shape functions gradients (to be re-implemented)
"""
function gradshfsps(reffe::LagrangianRefFE{D},points::Array{Point{D}}) where D
	grad  = Polynomials.gradient(reffe.prebasis,points)
	grad = reshape(grad,size(grad,1),size(grad,2)*size(grad,3))
	grad = reffe.changeofbasis*grad
	grad = reshape(grad,size(grad,1),:,D)
	nshf=size(reffe.changeofbasis,1)
	aux1=D^(reffe.rank+1)
	rgshfdim=[(i-1)*nshf+1:nshf*i for i=1:D^(reffe.rank)]
	gradshfss=zeros(Float64,rgshfdim[end][end],size(grad,2),aux1)
	for i=1:length(rgshfdim)
		for j=1:D
			k = D*(i-1)+j
			gradshfss[rgshfdim[i],:,k]=grad[:,:,j]
		end
	end
	aux2=D*ones(Int64,reffe.rank+1)
	gradshfss=reshape(gradshfss,size(gradshfss,1),size(gradshfss,2),aux2...)
end
