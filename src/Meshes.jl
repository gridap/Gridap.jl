# using ..Polytopes: PointInt

export gidscellxtype, gidscellxtypefast
export flip
export LexIndexSet
export lex2int
export cartesianindexmatrixoffset!
export Mesh, StructHexMesh

abstract type Mesh end

struct StructHexMesh <: Mesh
	# Here I am computing the cel vefs, and store the result in a big array, probably better to compute it on the fly on demand. Analogously for the dual mesh, for any relation of nface dim against nface dim. For the moment, I am creating a mesh using the same structures as unstructured meshes. A very light weight structured mesh in which arrays are replaced by functions, probably using an AbstractArray would be possible.
	cellvefs::Array{Array{Int64,1},1}
	vefcells::Array{Array{Int64,1},1}
	coordinates::Array{Float64,2}
	polytope::Polytope
	function StructHexMesh(nparts,nprocs=-1,iproc=0)
		D = length(nparts)
		if (iproc == 0)
			iproc = tuple(zeros(Int64,D)...)
		end
		extrusion = Polytopes.PointInt{D}(ones(Int64,D))
		polytope = Polytope(extrusion)
		numnfs=length(polytope.nfaces)
		nvefs=length(polytope.nfaces)-1
		aux = 2*ones(Int64,D)
		offst = [ prod(aux[1:i-1]) for i=1:D]
		# offst = linfs.offset
		nfspol = hcat([polytope.nfaces[i].extrusion for i=1:nvefs+1]...)
		# nftypep1 = lex2int(linfs,nfspol)
		nftype = [offst'*polytope.nfaces[i].extrusion for i=1:nvefs+1].+1
		spdimv=2*ones(Int64,D)
		exttypes = Array{Array{Int64,1}}(undef,spdimv...)
		cartesianindexmatrix!(exttypes)
		exttypes = hcat(reshape(exttypes, length(exttypes))...)'
		dimtypes = [sum(exttypes[i,:]) for i=1:size(exttypes,1)]
		extli = hcat([(nparts+flip.(exttypes[i,:])) for i=1:size(exttypes,1)]...)
		offnf= hcat([[prod(extli[1:i-1,j]) for i=1:D] for j=1:size(exttypes,1)]...)
		extli = [ prod(extli[1:D,j]) for j=1:size(exttypes,1)]
		perm = sortperm(dimtypes)
		invperm = sortperm(perm)
		offst = accumulate(+, extli[perm])
		offst[2:end]=offst[1:end-1]; offst[1]=0
		offst = offst[invperm]
		offstnf=offst[nftype]
		cells = Array{Array{Int64,1}}(undef,nparts...)
		cartesianindexmatrixoffset!(cells,iproc)
		cells=reshape(cells,length(cells))
		numcells=length(cells)
		cellvefs=gidscellxtype(cells,polytope,nparts)
		for ivef=1:nvefs+1
			cellvefs[ivef,:].+=offstnf[ivef]
		end
		cellvefs = [cellvefs[:,icell] for icell=1:size(cellvefs,2)]
		vefcells = [ Int64[] for i=1:offstnf[nvefs+1]+numcells]
		for cell=1:numcells
			for iface=1:nvefs+1
				push!(vefcells[cellvefs[cell][iface]],cell)
			end
		end
		# for cell=1:length(cells)
		# 	cellvefs[:,cell]+=offstnf
		# end
		coordinates = zeros(Float64,2,2)
		return new(cellvefs,vefcells,coordinates,polytope)
	end
end

struct LexIndexSet
    range::Vector{Int64}
    offset::Vector{Int64}
    function LexIndexSet(range::Vector{Int64})
        offset = [ prod(range[1:i-1]) for i=1:length(range)]
		return new(range,offset)
	end
end
function lex2int(liset::LexIndexSet,points::Array{Int64,1})
	return liset.offset'*points.+1
end

function gidscellxtype(cells,polytope,nparts)
	licells = LexIndexSet(nparts)
    numnfs=length(polytope.nfaces)
    gid1=zeros(Int64,numnfs)
    gid=zeros(Int64,numnfs,length(cells))
	offst = [nparts+flip.(polytope.nfaces[i].extrusion) for i=1:numnfs]
	for inface=1:numnfs
		offst[inface] = [ prod(offst[inface][1:i-1]) for i=1:length(offst[inface])]
		nface=polytope.nfaces[inface]
		anc=nface.anchor
		gid1[inface] = offst[inface]'*anc.+1
	end
	cellvefs = Array{Int64,2}(undef,numnfs,length(cells))
	for i=1:numnfs; cellvefs[i,:].=gid1[i]; end
	cellvefs+=hcat(offst...)'*hcat(cells...)
    return cellvefs
end

function flip(a)
    if (a==0)
        a=1
    elseif a==1
        a=0
    end
end

@generated function cartesianindexmatrixoffset!(A::Array{Array{Int64,M},N},offst) where {M,N}
    quote
		offst=offst.-1
        @nloops $N i A begin (@nref $N A i) = [map(+,(@ntuple $N i),offst)...] end
    end
end

@generated function cartesianindexmatrix!(A::Array{Array{T,M},N}) where {T,M,N}
  quote
    @nloops $N i A begin (@nref $N A i) = [((@ntuple $N i).-1)...] end
  end
end
