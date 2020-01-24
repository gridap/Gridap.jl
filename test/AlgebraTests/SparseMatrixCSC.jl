module SparseMatrixCSCTests

    using SparseArrays
    using LinearAlgebra
    using Gridap.Algebra
    using Test

    maxnz=5
    maxrows=5
    maxcols=5

    for T in (Int32,Int64,Float32,Float64)
        I = Vector{Int}()
        J = Vector{Int}()
        V = Vector{T}()
        for (ik, jk, vk) in zip(rand(1:maxrows, maxnz), rand(1:maxcols, maxnz), rand(1:T(maxnz), maxnz-1))
            push_coo!(I,J,V,ik,jk,vk)
        end
        push_coo!(I,J,V,maxrows,maxcols,maxnz)
        finalize_coo!(I,J,V,maxcols,maxrows)
        CSC = sparse(I, J, V, maxcols,maxrows)

        V2 = copy(V)
        V2 .= 0
        CSC2 = sparse(I, J, V2, maxcols,maxrows)

        copy_entries!(CSC2,CSC)
        @test CSC2 == CSC
        copy_entries!(CSC2,CSC2)
        @test CSC2 == CSC

        @test convert(SparseMatrixCSC, CSC) === CSC

        @test hasrowmajororder(CSC)  == false
        @test hascolmajororder(CSC) == true
        @test getptr(CSC)           == CSC.colptr
        @test getindices(CSC)       == rowvals(CSC)

        for col = 1:maxcols
            for j in nzrange(CSC, col)
                @test nonzeros(CSC)[j] == CSC[rowvals(CSC)[j], col]
            end
        end

        @test size(CSC) == (maxrows,maxcols)

        @test nnz(CSC) == count(!iszero, CSC) <= maxnz

        vold = getindex(CSC,maxrows,maxcols)
        add_entry!(CSC,1,maxrows,maxcols,+)
        @test getindex(CSC,maxrows,maxcols) == vold+1

        fill_entries!(CSC,0)
        @test all(x->x==0, nonzeros(CSC))



    end
    
end
