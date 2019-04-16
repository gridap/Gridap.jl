module FESpaces

export FESpace
export globalnumbering, computelgidvefs

using Numa.Helpers
using Numa.RefFEs
using Numa.Polytopes
using Numa.CellValues
using Numa.Geometry

using SparseArrays

using Numa.CellValues: CellVectorByComposition

# Abstract types and interfaces

"""
Abstract FE Space parameterized with respec to the environment dimension `D`,
the cell dimension `Z`, and the field type `T` (rank) and the number type `E`
"""
abstract type FESpace{D,Z,T,E} end

function globaldofs(::FESpace{D,Z,T,E}, cellvefs, vefcells)::CellVector{Int} where {D,Z,T,E}
  @abstractmethod
end

function get_grid(::FESpace{D,Z,T,E})::Grid{D,Z} where {D,Z,T,E}
	@abstractmethod
end

"""
Conforming FE Space, where only one RefFE is possible in the whole mesh
"""
struct ConformingFESpace{D,Z,T} <: FESpace{D,Z,T,Float64}
	# For the moment, I am not considering E (to think)
	reffe::LagrangianRefFE{D,T}
	grid::Grid{D,D}
end

# @santiagobadia : I would like to import grid... Check where is it being used
get_grid(this::ConformingFESpace) = this.grid

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
	grid = get_grid(this)
	graph = gridgraph(grid) # Generates the GridGraph associated with this grid.
	cellvefs = celltovefs(graph) # Show vefs for each cell
	vefcells = veftocells(graph) # Show cells around each vef
	# @santiagobadia : I had the celltovefs and veftocells calls inside globaldofs
	# Now here in order to call the previous methods once, since they can be
	# costly (for some concrete grids).
	gldofs = globaldofs(this, cellvefs, vefcells)
	ndofs = gldofs[end][end]
	# @santiagobadia : We are calling twice cell_to_vefs...
	cell_to_dofs = CellVectorByComposition(cellvefs, gldofs)
	ConformingAssembler{Float64}(cell_to_dofs, cell_to_dofs, ndofs)
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

function globaldofs(this::ConformingFESpace, cellvefs, vefcells)
	reffe = this.reffe
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

end # module FESpaces
