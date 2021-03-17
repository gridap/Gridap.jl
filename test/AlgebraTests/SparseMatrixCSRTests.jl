module SparseMatrixCSRTests

using SparseArrays
using SparseMatricesCSR
using LinearAlgebra
using Gridap.Algebra
using Test

maxnz=5
maxrows=5
maxcols=5
maxrowsorcols=7
int_types=(Int32,Int)
float_types=(Float32,Float64)
Bi_types=(0,1)

for Ti in int_types
  for Tv in float_types
    for Bi in Bi_types
      I = Vector{Ti}()
      J = Vector{Ti}()
      V = Vector{Tv}()
      for (ik, jk, vk) in zip(rand(1:maxrows, maxnz), rand(1:maxcols, maxnz), rand(1:Tv(maxnz), maxnz-1))
        push_coo!(SparseMatrixCSR{Bi,Tv,Ti},I,J,V,ik,jk,vk)
      end
      push_coo!(SparseMatrixCSR{Bi,Tv,Ti},I,J,V,maxrows,maxcols,maxnz)
      finalize_coo!(SparseMatrixCSR{Bi,Tv,Ti},I,J,V,maxrows,maxcols)
      CSC = sparse(I, J, V, maxrows,maxcols)
      CSR = sparsecsr(Val(Bi),I,J,V,maxrows,maxcols)

      @test is_entry_stored(SparseMatrixCSR{Bi,Tv,Ti},1,1)
      @test is_entry_stored(SparseMatrixCSR{Bi,Tv,Ti},1,2)
      @test is_entry_stored(SparseMatrixCSR{Bi,Tv,Ti},2,1)

      _l = 10
      _I, _J, _V = allocate_coo_vectors(SparseMatrixCSR{Bi,Tv,Ti},_l)
      @test length(_I) == _l
      @test length(_J) == _l
      @test length(_V) == _l
      @test eltype(_I) == Ti
      @test eltype(_J) == Ti
      @test eltype(_V) == Tv

      V2 = copy(V)
      V2 .= 0
      CSR2 = sparsecsr(Val(Bi),I, J, V2,maxrows,maxcols)
      copy_entries!(CSR2,CSR)
      @test CSR2 == CSR
      copy_entries!(CSR2,CSR2)
      @test CSR2 == CSR

      @test CSC == CSR
      @test nnz(CSC) == count(!iszero, CSC) == nnz(CSR) == count(!iszero, CSR)

      TCSC = sparse(J, I, V, maxrows, maxcols)
      TCSR = sparsecsr(Val(Bi), J, I, V, maxrows, maxcols)

      vold = getindex(CSR,maxrows,maxcols)
      add_entry!(+,CSR,1,maxrows,maxcols)
      @test getindex(CSR,maxrows,maxcols) == vold+1

      fill_entries!(CSR,0)
      @test all(x->x==0, nonzeros(CSR))

    end
  end
end

end # module
