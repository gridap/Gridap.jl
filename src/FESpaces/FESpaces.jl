module FESpaces

export FESpace
export ConformingFESpace
export globalnumbering, computelgidvefs
export interpolate
export Assembler
export ConformingAssembler

using Numa: evaluate
using Numa.Helpers
using Numa.RefFEs
using Numa.Polytopes
using Numa.CellValues
using Numa.Geometry

using SparseArrays

using Numa.CellValues: CellVectorByComposition
using Numa.CellMaps
using Numa.CellMaps: CellFieldFromExpand

import Numa: evaluate!

# Abstract types and interfaces

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank) and the number type `E`
"""
abstract type FESpace{D,Z,T,E} end

function assemblycellgids(::FESpace)::CellVector{Int}
	@abstractmethod
end

function applyconstraints(::FESpace,
	cellvec::CellVector)::Tuple{CellVector,CellVector{Int}}
	@abstractmethod
end

function applyconstraintsrows(::FESpace,
	cellmat::CellMatrix)::Tuple{CellMatrix,CellVector{Int}}
	@abstractmethod
end

function applyconstraintscols(::FESpace,
	cellmat::CellMatrix)::Tuple{CellMatrix,CellVector{Int}}
	@abstractmethod
end

"""
Conforming FE Space, where only one RefFE is possible in the whole mesh
"""
struct ConformingFESpace{D,Z,T} <: FESpace{D,Z,T,Float64}
	# For the moment, I am not considering E (to think)
	reffe::LagrangianRefFE{D,T}
	trian::Triangulation{D,Z}
	graph::GridGraph
	dof_eqclass::CellVector{Int}
	num_free_dofs::Int
	num_fixed_dofs::Int
	# dirichlet_tags::Array{Int} # Physical tags that describe strong dirichlet boundary
end

function ConformingFESpace(
	reffe::LagrangianRefFE{D,T},
	trian::Triangulation{D,Z},
	graph::GridGraph,
	# @santiagobadia : To be eliminated
	is_fixed_vef::AbstractVector{Bool}) where {D,Z,T}
	cellvefs = celltovefs(graph)
	vefcells = veftocells(graph)
	gldofs, nfree, nfixed  = globaldofs(reffe, cellvefs, vefcells, is_fixed_vef)
	# @santiagobadia : We are calling twice cell_to_vefs...
	dof_eqclass = CellVectorByComposition(cellvefs, gldofs)
	ConformingFESpace{D,Z,T}(reffe, trian, graph, dof_eqclass, nfree, nfixed)
end

function applyconstraints(this::ConformingFESpace,
	cellvec::CellVector)
	return cellvec, this.dof_eqclass
end

function applyconstraintsrows(this::ConformingFESpace,
	cellmat::CellMatrix)
	return cellmat, this.dof_eqclass
end

function applyconstraintscols(this::ConformingFESpace,
	cellmat::CellMatrix)
	return cellmat, this.dof_eqclass
end

"""
Abstract assembly operator
"""
abstract type Assembler{E} end

ndofs(::FESpace) = @abstractmethod

function assemblevector(::CellVector{E})::Vector{E} where {E}
	@abstractmethod
end

function assemblematrix(::CellMatrix{E})::Matrix{E} where {E}
	@abstractmethod
end

struct ConformingAssembler{E} <: Assembler{E}
	testfesp::FESpace
	trialfesp::FESpace
	num_dofs::Int
end

function ConformingAssembler(this::FESpace)
	ConformingAssembler{Int}(this, this, this.num_free_dofs)
end

function ConformingAssembler(test::FESpace, trial::FESpace)
	@assert trial.num_free_dofs == test.num_free_dofs
	ConformingAssembler{Int}(trial, test, trial.num_free_dofs)
end

# Methods

function assemble(this::Assembler, vals::CellVector{T}) where T
	_vals, rows_m = applyconstraints(this.testfesp, vals)
	# @santiagobadia : Evaluate efficiency, best way to do it in Julia
	# without pre-allocate loop?
	aux_row = []; aux_vals = []
	for (rows_c,vals_c) in zip(rows_m,_vals)
		for (i,gid) in enumerate(rows_c)
			if gid > 0
				aux_vals = [aux_vals..., vals_c[i]]
				aux_row = [ aux_row..., rows_c[i]...]
			end
		end
	end
	return Array(sparsevec(aux_row, aux_vals))
end

function assemble(this::Assembler, vals::CellMatrix{T}) where T
	_vals, rows_m = applyconstraintsrows(this.testfesp, vals)
	_vals, cols_m = applyconstraintscols(this.trialfesp, _vals, )
	# @santiagobadia : Evaluate efficiency, best way to do it in Julia
	# without pre-allocate loop?
	aux_row = []; aux_col = []; aux_vals = []
	for vals_c in _vals
	end
	for (rows_c, cols_c, vals_c) in zip(rows_m,cols_m,_vals)
		for (i,gidrow) in enumerate(rows_c)
			if gidrow > 0
				for (j,gidcol) in enumerate(cols_c)
					if gidcol > 0
						aux_row = [aux_row..., rows_c[i]]
						aux_col = [aux_col..., cols_c[j]]
						aux_vals = [aux_vals..., vals_c[i,j]]
						# for I in Iterators.product(rows_c, cols_c)
						# 	aux_row = [aux_row..., I[1]]
						# 	aux_col = [aux_col..., I[2]]
						# end
						# aux_vals = [aux_vals..., vec(vals_c)...]
					end
				end
			end
		end
	end
	return sparse(aux_row, aux_col, aux_vals)
end

function globaldofs(reffe::RefFE, cellvefs, vefcells, is_fixed_vef::AbstractVector)
	nfdofs=Array{Array{Int64},1}(undef,length(vefcells))
	c=1
	c_n = -1
	nfdofs_l = []
	nfdofs_g = zeros(Int, length(vefcells)+1)
	nfdofs_g[1] = 1
	for (ignf,nf) in enumerate(vefcells)
		owner_cell = nf[1]
		lid_vef = findfirst(i->i==ignf,cellvefs[owner_cell])
		num_nf_dofs = length(reffe.nfacedofs[lid_vef])
		if ( is_fixed_vef[ignf] )
			nfdofs_l = [nfdofs_l..., c_n:c_n-num_nf_dofs+1... ]
			c_n -= num_nf_dofs
		else
			nfdofs_l = [nfdofs_l..., c:c+num_nf_dofs-1... ]
			c += num_nf_dofs
		end
		nfdofs_g[ignf+1] += num_nf_dofs + nfdofs_g[ignf]
	end
	return CellVectorFromDataAndPtrs(nfdofs_l, nfdofs_g), c-1, -c_n-1
end

function interpolate(fun::Function, fesp::FESpace)
	reffe = fesp.reffe
	dofb = reffe.dofbasis
	trian = fesp.trian
	phi = geomap(trian)
	uphys = fun ∘ phi
	celldofs = fesp.dof_eqclass
	free_dofs = zeros(Float64, fesp.num_free_dofs)
	fixed_dofs = zeros(Float64, fesp.num_fixed_dofs)
	aux = zeros(Float64,cellsize(fesp.dof_eqclass)...)
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
	cdofs = CellVectorFromLocalToGlobalPosAndNeg(celldofs, free_dofs, fixed_dofs)
	# cdofs = CellVectorFromLocalToGlobal(celldofs,free_dofs)
	intu = CellFieldFromExpand(shb, cdofs)
end

end # module FESpaces
