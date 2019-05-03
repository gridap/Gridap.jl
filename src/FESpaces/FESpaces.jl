module FESpaces

export FESpace
export globalnumbering, computelgidvefs
export interpolate

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

# Abstract types and interfaces

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank) and the number type `E`
"""
abstract type FESpace{D,Z,T,E} end

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
end

function ConformingFESpace(
	reffe::LagrangianRefFE{D,T},
	trian::Triangulation{D,Z},
	graph::GridGraph) where {D,Z,T}
	cellvefs = celltovefs(graph)
	vefcells = veftocells(graph)
	gldofs = globaldofs(reffe, cellvefs, vefcells)
	ndofs = gldofs[end][end]
	# @santiagobadia : We are calling twice cell_to_vefs...
	dof_eqclass = CellVectorByComposition(cellvefs, gldofs)
	ConformingFESpace{D,Z,T}(reffe, trian, graph, dof_eqclass, ndofs, 0)
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
	assembly_op_rows::CellVector{Int}
	assembly_op_cols::CellVector{Int}
	num_dofs::Int
end

function ConformingAssembler(this::FESpace)
	ConformingAssembler{Int}(this.dof_eqclass, this.dof_eqclass, this.num_free_dofs)
end


# Methods

function assemble(this::Assembler, vals::CellVector{T}) where T
  rows_m = this.assembly_op_rows
	# @santiagobadia : Evaluate efficiency, best way to do it in Julia
	# without pre-allocate loop?
  aux_row = []; aux_vals = []
  for vals_c in vals
    aux_vals = [aux_vals..., vals_c...]
  end
  for rows_c in rows_m
    aux_row = [ aux_row..., rows_c...]
  end
  return Array(sparsevec(aux_row, aux_vals))
end

function assemble(this::Assembler, vals::CellMatrix{T}) where T
  rows_m = this.assembly_op_rows
  cols_m = this.assembly_op_cols
	# @santiagobadia : Evaluate efficiency, best way to do it in Julia
	# without pre-allocate loop?
  aux_row = []; aux_col = []; aux_vals = []
  for vals_c in vals
    aux_vals = [aux_vals..., vec(vals_c)...]
  end
  for (rows_c, cols_c) in zip(rows_m,cols_m)
    for I in Iterators.product(rows_c, cols_c)
      aux_row = [aux_row..., I[1]]
      aux_col = [aux_col..., I[2]]
    end
  end
  return sparse(aux_row, aux_col, aux_vals)
end

function globaldofs(reffe::RefFE, cellvefs, vefcells)
	nfdofs=Array{Array{Int64},1}(undef,length(vefcells))
	c=1
	nfdofs_l = []
	nfdofs_g = zeros(Int, length(vefcells)+1)
	nfdofs_g[1] = 1
	for (ignf,nf) in enumerate(vefcells)
		owner_cell = nf[1]
		lid_vef = findfirst(i->i==ignf,cellvefs[owner_cell])
		num_nf_dofs = length(reffe.nfacedofs[lid_vef])
		nfdofs_l = [nfdofs_l..., c:c+num_nf_dofs-1... ]
		nfdofs_g[ignf+1] += num_nf_dofs + nfdofs_g[ignf]
		c += num_nf_dofs
	end
	return CellVectorFromDataAndPtrs(nfdofs_l, nfdofs_g)
end

function interpolate(fun::Function, fesp::FESpace)
  reffe = fesp.reffe
  dofb = reffe.dofbasis
  trian = fesp.trian
  phi = geomap(trian)
  uphys = fun ∘ phi
  celldofs = fesp.dof_eqclass
  free_dofs = zeros(Float64, fesp.num_free_dofs)
  for (imap,l2g) in zip(uphys,celldofs)
    free_dofs[l2g] = evaluate(dofb,imap)
  end
  shb = ConstantCellValue(reffe.shfbasis, ncells(trian))
	cdofs = CellVectorFromLocalToGlobal(celldofs,free_dofs)
  intu = CellFieldFromExpand(shb, cdofs)
end

end # module FESpaces
