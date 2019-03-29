
export FESpace
export globalnumbering, computelgidvefs

using Numa.Meshes
using Numa.RefFEs
using Numa.Polytopes

"""
FE Space structure, where only one RefFE is possible in the whole mesh (to be improved in the future)
"""
struct FESpace
	reffe::LagrangianRefFE
	mesh::Mesh
	l2giddof::Array{Array{Int64,1},1}
	numgdof::Int64
	# Constructor not mutable HERE NOW
	function FESpace(reffe::RefFE,mesh::Mesh)
		giddof=globalnumbering(reffe,mesh)
		l2giddof=computelgidvefs(reffe,mesh,giddof[1])
		return new(reffe,mesh,l2giddof,giddof[2])
	end
end

function globalnumbering(reffe::RefFE,mesh::Mesh)
	nfdofs=Array{Array{Int64},1}(undef,length(mesh.vefcells))
	c=1
	igvef=0
	for igvef = 1:length(mesh.vefcells)
		# global c
		ownercell = mesh.vefcells[igvef][1]
		lidvef = findfirst(i->i==igvef,mesh.cellvefs[ownercell])
		numvefdofs =length(reffe.nodes.nfacenodes[lidvef])
		nfdofs[igvef] = [c:numvefdofs+c-1...]
		c+=numvefdofs
	end
	return [nfdofs,c-1]
end

function computelgidvefs(reffe::RefFE, mesh::Mesh,gldofs)
	polytope=reffe.polytope
	lgidvefs=Array{Array{Int64,1},1}(undef,length(mesh.cellvefs))
	for icell=1:length(mesh.cellvefs)
		lgidvefs[icell]= Array{Int64,1}(undef,size(reffe.nodes.coordinates,1))
		for lidvef = 1:length(polytope.nfaces)
			lgidvefs[icell][reffe.nodes.nfacenodes[lidvef]]=gldofs[mesh.cellvefs[icell][lidvef]]
		end
	end
	return lgidvefs
end

# Target:
# integrate(f[q]*inner(testspace[U].grad[q],trialspace[U].grad[q]))
# U = shapefunctions(V), V: fespace
# V = shapefunctions(W), W: fespace
# f = function
