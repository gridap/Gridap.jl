"""
Conforming FE Space, where only one RefFE is possible in the whole mesh
"""
struct ConformingFESpace{D,Z,T,E} <: FESpace{D,Z,T,E}
	# For the moment, I am not considering E (to think)
	reffes::LagrangianRefFE{D,T}
	triangulation::Triangulation{D,Z}
	nf_dofs::Vector{<:IndexCellArray{Int}}
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
	ConformingFESpace{D,Z,T,Float64}(reffe, trian, gldofs, cell_eqclass, nfree, nfixed, dir_tags)
end

for op in (:reffes, :triangulation, :nf_dofs, :cell_eqclass,
	:num_free_dofs, :num_fixed_dofs, :dir_tags)
	@eval begin
		$op(this::ConformingFESpace) = this.$op
	end
end

function applyconstraints(this::ConformingFESpace,
	cellvec::CellVector)
	return cellvec, this.nf_dofs
end

function applyconstraintsrows(this::ConformingFESpace,
	cellmat::CellMatrix)
	return cellmat, this.nf_dofs
end

function applyconstraintscols(this::ConformingFESpace,
	cellmat::CellMatrix)
	return cellmat, this.nf_dofs
end

ConformingFESpaces{D,Z,T,E} = Union{ConformingFESpace{D,Z,T,E}, FESpaceWithDirichletData{D,Z,T,E,ConformingFESpace{D,Z,T,E}}}
