export Polytope
export NodesArray
export NFace
export dim, numnftypes
#export nfaceperms, cartesianindexmatrix!, nodescoordinates, tensorfill!
#
"""
n-face of the polytope, i.e., any polytope of lower dimension (n) representing its boundary and the polytope itself (for n equal to the space dimension)
"""
mutable struct NFace
    anchor::Array{Int64,1}
    extrusion::Array{Int64,1}
end
dim(nface::NFace) = sum(nface.extrusion)

"""
Aggregation of all n-faces that compose the polytope boundary and the polytope itself, the classification of n-faces with respect to their dimension and type
"""
mutable struct Polytope
    extrusion::Array{Int64,1}
    nfaces::Array{NFace,1}
	dimnfs::Array{UnitRange{Int64},1}
	nftype::Array{Int64,1}
    function Polytope(extrusion::Array{Int64,1})
        polytope = new()
        polytope.extrusion = extrusion
        zeros = [0 for i=1:length(polytope.extrusion)]
        polytope.nfaces = []
        polytope.nfaces = createnfaces!(zeros, zeros, extrusion, polytope.nfaces)
        [sort!(polytope.nfaces, by = x -> x.anchor[i]) for i=1:length(polytope.extrusion)]
        [sort!(polytope.nfaces, by = x -> x.extrusion[i]) for i=1:length(polytope.extrusion)]
        [sort!(polytope.nfaces, by = x -> sum(x.extrusion))]
		# ) for i=1:length(polytope.nfaces)]
		numnfs = length(polytope.nfaces)
		nfsdim = [sum(polytope.nfaces[i].extrusion) for i=1:numnfs]
		spdims=length(polytope.extrusion)
		polytope.dimnfs = Array{UnitRange{Int64},1}(undef,spdims+1)
		dim=0; i=1
		for iface=1:numnfs
			if (nfsdim[iface]>dim)
				# global dim
				# global i
				dim+=1
				polytope.dimnfs[dim]=i:iface-1
				i=iface
			end
		end
		polytope.dimnfs[spdims+1]=numnfs:numnfs
		aux=2*ones(Int64,spdims)
		offst=[ prod(aux[1:i-1]) for i=1:spdims]
		polytope.nftype=[offst'*polytope.nfaces[i].extrusion for i=1:length(polytope.nfaces)].+1
        return polytope
    end
end

dim(polytope::Polytope) = length(polytope.extrusion)
numnftypes(polytope::Polytope) = 2^dim(polytope)
# numnfs(polytope::Polytope) = length(polytope.nfaces)
# numvefs(polytope::Polytope) = length(polytope.nfaces)-1


# Create list of all nfaces of a polytope sorted by dim first
function createnfaces!(anchor::Array{Int64,1}, extrusion::Array{Int64,1}, extend::Array{Int64,1}, list)
    newext = copy(extend) # extend cannot be modified for recursion
    numdims = length(newext)
    a = NFace(anchor, extrusion)
    list = [list...,a]
    for i = 1:numdims
        if (newext[i] > 0)
            edim = [0 for j=1:numdims]; edim[i] = 1; newext[i]=0
            list = createnfaces!(anchor+edim, extrusion, newext, list)
            list = createnfaces!(anchor, extrusion+edim, newext, list)
        end
    end
    return list
end

"""
Array of nodes for a give polytope and order
"""
mutable struct NodesArray
    coordinates::Array{Float64,N} where N
    nfacenodes::Array{Array{Int64,1},1}
    closurenfacenodes::Array{Array{Int64,1},1}
    function NodesArray(polytope::Polytope,orders::Array{Int64,1})
        nodes=new()
        nodes.closurenfacenodes = [createnodes(polytope.nfaces[i], orders) for i=1:length(polytope.nfaces)]
        nodes.nfacenodes = [createnodes(polytope.nfaces[i], orders, isopen=true) for i=1:length(polytope.nfaces)]
        spdims = length(orders)
        coords = [nodescoordinates(orders[i], nodestype="Equispaced") for i=1:spdims]
        ordsp1=orders.+1; dims = tuple(ordsp1...)
        nodes.coordinates=Array{Float64,spdims+1}(undef,tuple(dims...,spdims))
        tensorfill!(nodes.coordinates,coords)
        nodes.coordinates = reshape(nodes.coordinates,:,spdims)
        return nodes
    end
end

# Create list of nface nodes with polytope indexing
function createnodes(nface::NFace,orders; isopen::Bool=false)
    spdims=length(nface.extrusion)
    @assert spdims == length(orders) "nface and orders dim must be identical"
    perms=nfaceperms(nface)
    p=perms[1]; pinv=perms[2]
    ordnf=orders[pinv]; ordp1=ordnf.+((isopen) ? d=-1 : d=1)
    A = Array{Array{Int64,1}}(undef,ordp1...)
    cartesianindexmatrix!(A)
    A = hcat(reshape(A, length(A))...)'
    if (isopen) A.+=1 end
    B = Array{Int64,2}(undef,prod(ordp1),spdims)
    for i=1:spdims
        if (p[i]==0) B[:,i].=nface.anchor[i]*orders[i] end
    end
    B[:,pinv] = copy(A)
    offst = orders.+1; offst = [ prod(offst[1:i-1]) for i=1:length(offst)]
    return B*offst.+1
end

# nface to polytope cartesian index change of basis and inverse
function nfaceperms(nface)
    pd = length(nface.extrusion)
    e = nface.extrusion
    c=0; p=[]
    for i=1:pd
        (e[i]!=0) ? (c+=1; p=[p...,c]) : p=[p...,0]
    end
    c=0; pinv=[]
    for i=1:pd
        (p[i]!=0) ? (c+=1; pinv=[pinv...,i]) : 0
    end
    return [p,pinv]
end

# 1-dim node coordinates
function nodescoordinates(order::Int; nodestype::String="Equispaced")
    @assert ((nodestype=="Equispaced") | (nodestype=="Chebyshev")) "Node type not implemented"
    ordp1 = order+1
    if nodestype == "Chebyshev"
        nodescoordinates = [ cos((i-1)*pi/order) for i=1:ordp1]
    elseif nodestype == "Equispaced"
        (order != 0) ? nodescoordinates = (2/order)*[ i-1 for i=1:ordp1].-1 : nodescoordinates = [0.0]
    end
    return nodescoordinates
end

@generated function cartesianindexmatrix!(A::Array{Array{T,M},N}) where {T,M,N}
    quote
        @nloops $N i A begin (@nref $N A i) = [((@ntuple $N i).-1)...] end
    end
end

@generated function tensorfill!(A::Array{T,N},c) where {T,N}
    quote
        @nloops $N i A begin
            t = @ntuple $N i; d = length(t)
            dim = t[d]
            nodedim = t[dim]
            (@nref $N A i) = c[dim][nodedim]
        end
    end
end

# Example of recursive type
#export NFace
#Export MyPolytope
#
#Mutable struct MyPolytope{R<:Integer}
#    dim :: R
#    nfaces :: MyPolytope{R}
#    function MyPolytope(dim::R) where R<:Integer
#      time = new{R}()
#      time.dim = dim
#      time
#    end
#End
#
#Function createnfaces(dim::R) where R<:Integer
#  a = MyPolytope(dim)
#  if (dim > 0)
#    a.nfaces = createnfaces(dim-1)
#  else
#    a.nfaces = a
#  end
#  return a
#End
