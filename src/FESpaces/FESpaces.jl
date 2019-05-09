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

reffes(::FESpace) = @abstractmethod

triangulation(::FESpace) = @abstractmethod

nf_eqclass(::FESpace) = @abstractmethod

cell_eqclass(::FESpace) = @abstractmethod

num_free_dofs(::FESpace) = @abstractmethod

num_fixed_dofs(::FESpace) = @abstractmethod

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
	reffes::LagrangianRefFE{D,T}
	triangulation::Triangulation{D,Z}
	nf_eqclass::Vector{<:IndexCellArray{Int}}
	cell_eqclass::IndexCellArray{Int}
	num_free_dofs::Int
	num_fixed_dofs::Int
	dir_tags::NTuple{M,Int} where M
end

function ConformingFESpace(
	reffe::LagrangianRefFE{D,T},
	trian::Triangulation{D,Z},
	graph::FullGridGraph,
	labels::FaceLabels) where {D,Z,T}
	return ConformingFESpace(reffe, trian, graph, labels, ())
end

function ConformingFESpace(
	reffe::LagrangianRefFE{D,T},
	trian::Triangulation{D,Z},
	graph::FullGridGraph,
	labels::FaceLabels,
	dir_tags::NTuple{M,Int}) where {D,Z,T,M}
	gldofs, nfree, nfixed  = globaldofs(reffe, graph, labels, dir_tags)
	cellvefs_dim = [connections(graph,D,i) for i in 0:D]
	offset = length.(gldofs)
	for i in 2:length(offset)
		offset[i] += offset[i-1]
	end
	offset = tuple(offset...)
	cellvefs = IndexCellValueByLocalAppendWithOffset(offset, cellvefs_dim...)
	dofs_all = IndexCellValueByGlobalAppend(gldofs...)
	cell_eqclass = CellVectorByComposition(cellvefs, dofs_all)
	ConformingFESpace{D,Z,T}(reffe, trian, gldofs, cell_eqclass, nfree, nfixed, dir_tags)
end

for op in (:reffes, :triangulation, :nf_eqclass, :cell_eqclass,
	:num_free_dofs, :num_fixed_dofs, :dir_tags)
	@eval begin
		$op(this::ConformingFESpace) = this.$op
	end
end

function applyconstraints(this::ConformingFESpace,
	cellvec::CellVector)
	return cellvec, this.nf_eqclass
end

function applyconstraintsrows(this::ConformingFESpace,
	cellmat::CellMatrix)
	return cellmat, this.nf_eqclass
end

function applyconstraintscols(this::ConformingFESpace,
	cellmat::CellMatrix)
	return cellmat, this.nf_eqclass
end

struct FESpaceWithDirichletData{D,Z,T,E,V<:FESpace{D,Z,T,E}} <: FESpace{D,Z,T,E}
	fesp::V
	dir_data::Vector{Float64}
end

for op in (:reffes, :triangulation, :gridgraph, :nf_eqclass, :cell_eqclass,
	:num_free_dofs, :num_fixed_dofs, :dir_tags)
	@eval begin
		$op(this::FESpaceWithDirichletData) = $op(this.fesp)
	end
end

dir_data(this::FESpaceWithDirichletData) = this.dir_data

function TestFESpace(this::FESpace)
  dv = zeros(Float64,num_fixed_dofs(this))
  return FESpaceWithDirichletData(this, dv)
end

function TrialFESpace( this::FESpace{D}, fun::Vector{Function}, labels::FaceLabels) where {D}
	dv = interpolate_dirichlet_data(fun, this, labels)
	# @santiagobadia : Put labels in FESPace
  return FESpaceWithDirichletData(this, dv)
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
	ConformingAssembler{Int}(this, this, num_free_dofs(this))
end

function ConformingAssembler(test::FESpace, trial::FESpace)
	@assert trial.num_free_dofs == test.num_free_dofs
	ConformingAssembler{Int}(trial, test, num_free_dofs(trial))
end

# Methods

function assemble(this::Assembler, vals::CellVector{T}) where T
	_vals, rows_m = applyconstraints(this.testfesp, vals)
	# @santiagobadia : Evaluate efficiency, best way to do it in Julia
	# without pre-allocate loop?
	aux_row = Int[]; aux_vals = Int[]
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
	aux_row = Int[]; aux_col = Int[]; aux_vals = Int[]
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
					end
				end
			end
		end
	end
	return sparse(aux_row, aux_col, aux_vals)
end

function globaldofs(reffe::RefFE{D,T},
	gridgr::FullGridGraph,
	labels::FaceLabels,
	dirt::NTuple{N,Int}) where {D,T,N}
	dim_eqclass = Int[]
	c=1
	c_n = -1
	for vef_dim in 0:D
		vefcells= connections(gridgr,vef_dim,D)
		cellvefs= connections(gridgr,D,vef_dim)
		vef_labels = labels_on_dim(labels,vef_dim)
		num_vefs = length(vefcells)
		nfdofs=Array{Array{Int64},1}(undef,num_vefs)
		nfdofs_l = Int[]
		nfdofs_g = zeros(Int, num_vefs+1)
		nfdofs_g[1] = 1
		for (ignf,nf) in enumerate(vefcells)
			owner_cell = nf[1]
			lid_vef_dim = findfirst(i->i==ignf,cellvefs[owner_cell])
			lid_vef = reffe.polytope.nf_dim[end][vef_dim+1][1]+lid_vef_dim-1
			# @santiagobadia : Better a method for nfs of a particular type...
			num_nf_dofs = length(reffe.nfacedofs[lid_vef])
			if ( _is_fixed(vef_labels[ignf],dirt,labels) )
				nfdofs_l = Int[nfdofs_l..., c_n:-1:c_n-num_nf_dofs+1... ]
				c_n -= num_nf_dofs
			else
				nfdofs_l = Int[nfdofs_l..., c:c+num_nf_dofs-1... ]
				c += num_nf_dofs
			end
			nfdofs_g[ignf+1] += num_nf_dofs + nfdofs_g[ignf]
		end
		dim_eqclass = [ dim_eqclass..., CellVectorFromDataAndPtrs(nfdofs_l, nfdofs_g) ]
	end
	return [ dim_eqclass , c-1, -c_n-1 ]
end

function interpolate(fun::Function, fesp::FESpace{D}) where {D}
	reffe = reffes(fesp)
	dofb = reffe.dofbasis
	trian = triangulation(fesp)
	phi = geomap(trian)
	uphys = fun ∘ phi
	celldofs = cell_eqclass(fesp)
	nfdofs = nf_eqclass(fesp)
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
end

function interpolate_dirichlet_data(fun::Vector{Function}, fesp::FESpace{D}, labels::FaceLabels) where {D}
	nf_labs_all = [ labels_on_dim(labels,idim) for idim in 0:D]
	nf_dofs_all = nf_eqclass(fesp)
	dtags = dir_tags(fesp)
	fixed_dofs = zeros(Float64, num_fixed_dofs(fesp))
	for (ifunc,f) in enumerate(fun)
		fh = interpolate(f, fesp)
		# Implement a new interpolate restricted to cells on the boundary for performance
		for idim in 0:D
			nf_labs = nf_labs_all[idim+1]
			nf_dofs = nf_dofs_all[idim+1]
			fh_fixed_dofs = fh.coeffs.gid_to_val_neg
			# How to extract this part? Do it correctly, with a FEFunction
			for (nf,nflab) in enumerate(nf_labs)
				if (_is_fixed(nflab, (dtags[ifunc],), labels))
					for dof in nf_dofs[nf]
						dof *= -1
						fixed_dofs[dof] = fh_fixed_dofs[dof]
					end
				end
			end
		end
	end
	return fixed_dofs
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

end # module FESpaces
