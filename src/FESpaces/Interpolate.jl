function interpolate(fun::Function, fesp::ConformingFESpaces{D}) where {D}
# function interpolate(fun::Function, fesp::FESpace{D}, ::ConformingMesh) where {D}
	reffe = _reffes(fesp)
	dofb = reffe.dofbasis
	trian = _triangulation(fesp)
	phi = geomap(trian)
	uphys = fun ∘ phi
	celldofs = cell_eqclass(fesp)
	nfdofs = nf_dofs(fesp)
	maxs = max([length(nfdofs[i]) for i=1:D+1]...)
	free_dofs = zeros(Float64, num_free_dofs(fesp))
	fixed_dofs = zeros(Float64, num_fixed_dofs(fesp))
	aux = zeros(Float64, maxs)
	for (imap,l2g) in zip(uphys,celldofs)
		evaluate!(dofb,imap,aux)
		for (i,gdof) in enumerate(l2g)
			if (gdof > 0)
				free_dofs[gdof] = aux[i]
			else
				fixed_dofs[-gdof] = aux[i]
			end
		end
	end
	shb = ConstantCellValue(reffe.shfbasis, ncells(trian))
	if (typeof(fesp) <: FESpaceWithDirichletData)
		cdofs = CellVectorFromLocalToGlobalPosAndNeg(celldofs, free_dofs, dir_data(fesp) )
	else
		cdofs = CellVectorFromLocalToGlobalPosAndNeg(celldofs, free_dofs, fixed_dofs)
	end
	intu = CellFieldFromExpand(shb, cdofs)
	return ConformingFEFunction(fesp, intu)
end

function interpolate_dirichlet_data(fun::Vector{Function}, fesp::FESpace{D}) where {D}
	labels = _labels(fesp)
	nf_labs_all = [ labels_on_dim(labels,idim) for idim in 0:D]
	nf_dofs_all = nf_dofs(fesp)
	dtags = dir_tags(fesp)
	@assert length(dtags) == length(fun)
	fixed_dofs_all = zeros(Float64, num_fixed_dofs(fesp))
	for (ifunc,f) in enumerate(fun)
		fh = interpolate(f, fesp)
		# Implement a new interpolate restricted to cells on the boundary for performance
		for idim in 0:D
			nf_labs = nf_labs_all[idim+1]
			nf_dofs = nf_dofs_all[idim+1]
			fh_fixed_dofs = fixed_dofs(fh)
			# How to extract this part? Do it correctly, with a FEFunction
			for (nf,nflab) in enumerate(nf_labs)
				if (_is_fixed(nflab, (dtags[ifunc],), labels))
					for dof in nf_dofs[nf]
						dof *= -1
						fixed_dofs_all[dof] = fh_fixed_dofs[dof]
					end
				end
			end
		end
	end
	return fixed_dofs_all
end

@inline function _is_fixed(v,dt,labels)
	for tag in dt
		labs = labels_on_tag(labels,tag)
		for label in labs
			if (label == v)
				return true
			end
		end
	end
	return false
end
