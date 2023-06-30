module BlockSparseMatrixAssemblersTests
using Test, BlockArrays, SparseArrays, LinearAlgebra

using Gridap
using Gridap.FESpaces, Gridap.ReferenceFEs, Gridap.MultiField

import Gridap.Algebra: nz_counter, nz_allocation, create_from_nz
import Gridap.Algebra: ArrayBuilder, SparseMatrixBuilder
import Gridap.Arrays: TouchEntriesMap, AddEntriesMap
import Gridap.Fields: MatrixBlock, VectorBlock, ArrayBlock

sol(x) = sum(x)

model = CartesianDiscreteModel((0.0,1.0,0.0,1.0),(5,5))
Ω = Triangulation(model)

reffe = LagrangianRefFE(Float64,QUAD,1)
V = FESpace(Ω, reffe; dirichlet_tags="boundary")
U = TrialFESpace(sol,V)

dΩ = Measure(Ω, 2)
biform((u1,u2),(v1,v2)) = ∫(∇(u1)⋅∇(v1) + u2⋅v2 - u1⋅v2)*dΩ
liform((v1,v2)) = ∫(v1 - v2)*dΩ

############################################################################################
# Normal assembly 

Y = MultiFieldFESpace([V,V])
X = MultiFieldFESpace([U,U])

u = get_trial_fe_basis(X)
v = get_fe_basis(Y)

data = collect_cell_matrix_and_vector(X,Y,biform(u,v),liform(v))
matdata = collect_cell_matrix(X,Y,biform(u,v))
vecdata = collect_cell_vector(Y,liform(v))  

assem = SparseMatrixAssembler(X,Y)
mat_builder = get_matrix_builder(assem)

m1 = nz_counter(get_matrix_builder(assem),(get_rows(assem),get_cols(assem)))
symbolic_loop_matrix!(m1,assem,matdata)
m2 = nz_allocation(m1)
numeric_loop_matrix!(m2,assem,matdata)
m3 = create_from_nz(m2)


v1 = nz_counter(get_vector_builder(assem),(get_rows(assem),))
symbolic_loop_vector!(v1,assem,vecdata)
v2 = nz_allocation(v1)
numeric_loop_vector!(v2,assem,vecdata)
v3 = create_from_nz(v2)


############################################################################################
# Block MultiFieldStyle

Gridap.Algebra.LoopStyle(a::ArrayBlock) = Gridap.Algebra.LoopStyle(first(a.array))

function Gridap.Algebra.nz_counter(builder::MatrixBlock,axis)
  s = size(builder.array)
  rows = map(i->axis[1][Block(i)],1:s[1])
  cols = map(i->axis[2][Block(i)],1:s[2])
  counters = [nz_counter(builder.array[i,j],(rows[i],cols[j])) for i in 1:s[1], j in 1:s[2]]
  return Gridap.Fields.ArrayBlock(counters,fill(true,size(counters)))
end

function Gridap.Algebra.nz_counter(builder::VectorBlock,axis)
  s = size(builder.array)
  rows = map(i->axis[1][Block(i)],1:s[1])
  counters = [nz_counter(builder.array[i],(rows[i],)) for i in 1:s[1]]
  return Gridap.Fields.ArrayBlock(counters,fill(true,size(counters)))
end

function Gridap.Algebra.nz_allocation(a::ArrayBlock)
  array = map(Gridap.Algebra.nz_allocation,a.array)
  return ArrayBlock(array,a.touched)
end

function Gridap.Algebra.create_from_nz(a::ArrayBlock)
  println("create_from_nz")
  array = map(Gridap.Algebra.create_from_nz,a.array)
  return mortar(array)
end

for T in (:AddEntriesMap,:TouchEntriesMap)
  @eval begin

    function Gridap.Fields.return_cache(
      k::$T,A::MatrixBlock,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
      qs = findall(v.touched)
      i, j = Tuple(first(qs))
      cij = return_cache(k,A.array[i,j],v.array[i,j],I.array[i],J.array[j])
      ni,nj = size(v.touched)
      cache = Matrix{typeof(cij)}(undef,ni,nj)
      for j in 1:nj
        for i in 1:ni
          if v.touched[i,j]
            cache[i,j] = return_cache(k,A.array[i,j],v.array[i,j],I.array[i],J.array[j])
          end
        end
      end
      cache
    end

    function Gridap.Fields.evaluate!(
      cache, k::$T,A::MatrixBlock,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
      ni,nj = size(v.touched)
      for j in 1:nj
        for i in 1:ni
          if v.touched[i,j]
            evaluate!(cache[i,j],k,A.array[i,j],v.array[i,j],I.array[i],J.array[j])
          end
        end
      end
    end

    function Gridap.Fields.return_cache(
      k::$T,A::VectorBlock,v::VectorBlock,I::VectorBlock)

      qs = findall(v.touched)
      i = first(qs)
      ci = return_cache(k,A.array[i],v.array[i],I.array[i])
      ni = length(v.touched)
      cache = Vector{typeof(ci)}(undef,ni)
      for i in 1:ni
        if v.touched[i]
          cache[i] = return_cache(k,A.array[i],v.array[i],I.array[i])
        end
      end
      cache
    end

    function Gridap.Fields.evaluate!(
      cache, k::$T,A::VectorBlock,v::VectorBlock,I::VectorBlock)
      ni = length(v.touched)
      for i in 1:ni
        if v.touched[i]
          evaluate!(cache[i],k,A.array[i],v.array[i],I.array[i])
        end
      end
    end

  end
end

mfs = BlockMultiFieldStyle()
Yb = MultiFieldFESpace([V,V];style=mfs)
Xb = MultiFieldFESpace([U,U];style=mfs)

ub = get_trial_fe_basis(Xb)
vb = get_fe_basis(Yb)

bdata = collect_cell_matrix_and_vector(Xb,Yb,biform(ub,vb),liform(vb))
bmatdata = collect_cell_matrix(Xb,Yb,biform(ub,vb))
bvecdata = collect_cell_vector(Yb,liform(vb))

block_assem = SparseMatrixAssembler(Xb,Yb)
block_assemblers = block_assem.block_assemblers
block_matrix_builders = ArrayBlock(map(get_matrix_builder,block_assemblers),fill(true,size(block_assemblers)))
block_vector_builders = ArrayBlock(map(get_vector_builder,diag(block_assemblers)),fill(true,2))
dummy_assem = block_assemblers[1]

bm1 = nz_counter(block_matrix_builders,(get_rows(block_assem),get_cols(block_assem)))
symbolic_loop_matrix!(bm1,dummy_assem,bmatdata)
bm2 = nz_allocation(bm1)
numeric_loop_matrix!(bm2,dummy_assem,bmatdata)
bm3 = create_from_nz(bm2)

bv1 = nz_counter(block_vector_builders,(get_rows(block_assem),))
symbolic_loop_vector!(bv1,dummy_assem,bvecdata)
bv2 = nz_allocation(bv1)
numeric_loop_vector!(bv2,dummy_assem,bvecdata)
bv3 = create_from_nz(bv2)

bv3 ≈ v3
bm3 ≈ m3

end # module