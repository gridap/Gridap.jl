export LagrangianRefFE
export shfsps, gradshfsps

"""
Abstract Reference Finite Element
"""
abstract type RefFE end

"""
Reference Finite Element a la Ciarlet, i.e., it relies on a local function (polynomial) space, an array of nodes (DOFs), and a polytope (cell topology). The rank of the approximating field can be arbitrary. The current implementation relies on the prebasis (e.g., monomial basis of polynomials) and a change-of-basis (using the node array) to generate the canonical basis, i.e., the shape functions.
"""
struct LagrangianRefFE<:RefFE
	polytope::Polytope
	prebasis::TensorProductPolynomialBasis
    nodes::NodesArray
	changeofbasis::Array{Float64,2}
	#dofs::Array{Function,1} (Other way to express it... as a function?)
	# fieldrank:: scalar, vector, tensor (or better field rank)
	# conformity, continuity (here or in the DOF numbering algorithm?)
	# dofsnface, dofsclosurenface
	dofsnface::Array{Array{Int64},1}
	rank::Int64
	numdof::Int64
    function LagrangianRefFE(polytope::Polytope,orders::Array{Int64,1},rank::Int64=0)
		nodes=NodesArray(polytope,orders)
		prebasis=TensorProductPolynomialBasis(orders, basistype="Monomial")
		dofspb=prebasis(nodes.coordinates)
		# dofspbs = dofbasis(prebasis)
		changeofbasis=inv(dofspb)
		dofsnface=nodes.nfacenodes
		spdim=dim(polytope)
		numdof = size(changeofbasis,1)*(spdim)^rank
		reffe=new(polytope, prebasis, nodes, changeofbasis, dofsnface, rank, numdof)
        return reffe
    end
end

"""
Compute shape functions in a set of points
"""
function shfsps(reffe::LagrangianRefFE,points)
	spdim=dim(reffe.polytope)
	scalm=reffe.changeofbasis*reffe.prebasis(points)
	nshf=size(reffe.changeofbasis,1)
	rgshfdim=[(i-1)*nshf+1:nshf*i for i=1:spdim^(reffe.rank)]
	aux1=dim(reffe.polytope)^reffe.rank
	shfss=zeros(Float64,rgshfdim[end][end],size(scalm,2),aux1)
	for i=1:length(rgshfdim)
		shfss[rgshfdim[i],:,i]=scalm
	end
	aux2=dim(reffe.polytope)*ones(Int64,reffe.rank)
	shfss=reshape(shfss,size(shfss,1),size(shfss,2),aux2...)
	return shfss
end

"""
Compute shape functions gradients (to be re-implemented)
"""
function gradshfsps(reffe::LagrangianRefFE,points)
	spdim = length(reffe.prebasis.polynomials)
	grad  = gradient(reffe.prebasis,points)
	grad = reshape(grad,size(grad,1),size(grad,2)*size(grad,3))
	grad = reffe.changeofbasis*grad
	grad = reshape(grad,size(grad,1),:,spdim)
	nshf=size(reffe.changeofbasis,1)
	aux1=dim(reffe.polytope)^(reffe.rank+1)
	rgshfdim=[(i-1)*nshf+1:nshf*i for i=1:spdim^(reffe.rank)]
	gradshfss=zeros(Float64,rgshfdim[end][end],size(grad,2),aux1)
	for i=1:length(rgshfdim)
		for j=1:spdim
			k = spdim*(i-1)+j
			gradshfss[rgshfdim[i],:,k]=grad[:,:,j]
		end
	end
	aux2=dim(reffe.polytope)*ones(Int64,reffe.rank+1)
	gradshfss=reshape(gradshfss,size(gradshfss,1),size(gradshfss,2),aux2...)
end
