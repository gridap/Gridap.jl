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
