function globaldofs(this::ConformingFESpace)
	mesh = this.mesh
	reffe = this.reffe
	nfdofs=Array{Array{Int64},1}(undef,length(mesh.vefcells))
	c=1
	nfdofs_l = []
	nfdofs_g = zeros(Int, length(mesh.vefcells)+1)
	nfdofs_g[1] = 1
	for (ignf,nf) in enumerate(mesh.vefcells)
		owner_cell = nf[1]
		lid_vef = findfirst(i->i==ignf,mesh.cellvefs[owner_cell])
		num_nf_dofs = length(reffe.nfacedofs[lid_vef])
		nfdofs_l = [nfdofs_l..., c:c+num_nf_dofs-1... ]
		nfdofs_g[ignf+1] += num_nf_dofs + nfdofs_g[ignf]
		c += num_nf_dofs
	end
	nfdofs_g
	nfdofs_l
	return CellVectorFromDataAndPtrs(nfdofs_l, nfdofs_g)
end

# Not sure we need this part... we could work vef based... and create a cell
# array from the vefarray
# function computelgidvefs(reffe::RefFE, mesh::Mesh,gldofs)
# 	polytope=reffe.polytope
# 	lgidvefs=Array{Array{Int64,1},1}(undef,length(mesh.cellvefs))
# 	for icell=1:length(mesh.cellvefs)
# 		lgidvefs[icell]= Array{Int64,1}(undef,size(reffe.nodes.coordinates,1))
# 		for lidvef = 1:length(polytope.nfaces)
# 			lgidvefs[icell][reffe.nodes.nfacenodes[lidvef]]=gldofs[mesh.cellvefs[icell][lidvef]]
# 		end
# 	end
# 	return lgidvefs
# end
