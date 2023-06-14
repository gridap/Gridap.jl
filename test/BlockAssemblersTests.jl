using Gridap
using Gridap.FESpaces, Gridap.Geometry, Gridap.CellData, Gridap.ReferenceFEs, Gridap.Fields

import Gridap.FESpaces: nz_counter, nz_allocation, create_from_nz
import Gridap.FESpaces: map_cell_cols, map_cell_rows

using BlockArrays, SparseArrays

############################################################################################

############################################################################################

sol(x) = sum(x)

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(5,5))
Ω = Triangulation(model)

reffe = LagrangianRefFE(Float64,QUAD,1)
V = FESpace(Ω, reffe; dirichlet_tags="boundary")
U = TrialFESpace(sol,V)

Y = MultiFieldFESpace([V,V])
X = MultiFieldFESpace([U,U])

dΩ = Measure(Ω, 2)
biform((u1,u2),(v1,v2)) = ∫(∇(u1)⋅∇(v1) + u2⋅v2)*dΩ
liform((v1,v2)) = ∫(v1 - v2)*dΩ

op = AffineFEOperator(biform,liform,X,Y)

u = get_trial_fe_basis(X)
v = get_fe_basis(Y)
data = collect_cell_matrix_and_vector(X,Y,biform(u,v),liform(v))
matdata = collect_cell_matrix(X,Y,biform(u,v))


struct BlockSparseMatrixAssembler <: Gridap.FESpaces.SparseMatrixAssembler
  glob_assembler   :: SparseMatrixAssembler
  block_assemblers :: AbstractArray{<:SparseMatrixAssembler}
  function BlockSparseMatrixAssembler(X::MultiFieldFESpace,Y::MultiFieldFESpace)
    nblocks = length(X)
    Gridap.Helpers.@check nblocks == length(Y)
    glob_assembler   = SparseMatrixAssembler(X,Y)

    row_offsets = [0,get_rows(glob_assembler).lasts...]
    col_offsets = [0,get_cols(glob_assembler).lasts...]

    block_assemblers = Matrix{SparseMatrixAssembler}(undef,nblocks,nblocks)
    for i in 1:nblocks
      for j in 1:nblocks
        row_map(row) = row .- row_offsets[i]
        col_map(col) = col .- col_offsets[j]
        row_mask(row) = true
        col_mask(col) = true
        strategy = GenericAssemblyStrategy(row_map,col_map,row_mask,col_mask)
        T = get_dof_value_type(X[i])
        matrix_type = SparseMatrixCSC{T,Int}
        vector_type = Vector{T}
        block_assemblers[i,j] = SparseMatrixAssembler(matrix_type,vector_type,X[i],Y[j],strategy)
      end
    end
    new{}(glob_assembler,block_assemblers)
  end
end

for fun in [:get_rows,:get_cols,:get_matrix_builder,:get_vector_builder,:get_assembly_strategy]
  @eval begin
    function Gridap.FESpaces.$fun(a::BlockSparseMatrixAssembler)
      $fun(a.glob_assembler)
    end
  end
end

"""
  TODO: We need to detect inactive blocks and avoid assembling them.
  Otherwise, we allocate unnecessary memory.
"""
function Gridap.FESpaces.assemble_matrix(ba::BlockSparseMatrixAssembler,matdata)

  rows = get_rows(ba.glob_assembler)
  cols = get_cols(ba.glob_assembler)
  r = rows.lasts .- [0,rows.lasts[1:end-1]...]
  c = cols.lasts .- [0,cols.lasts[1:end-1]...]
  A = BlockMatrix{Float64}(undef_blocks,r,c)

  block_assemblers = ba.block_assemblers
  for i in 1:length(r)
    for j in 1:length(c)
      a = block_assemblers[i,j]
      _matdata = (map(y->lazy_map(x->getindex(x,i,j),y),matdata[1]),
                  map(y->lazy_map(x->getindex(x,i),y),matdata[2]),
                  map(y->lazy_map(x->getindex(x,j),y),matdata[3])
                  )
      A[Block(i,j)] = assemble_matrix(a,_matdata)
    end
  end

  return A
end


ba = BlockSparseMatrixAssembler(X,Y)
mat_blocks = assemble_matrix(ba,matdata)
