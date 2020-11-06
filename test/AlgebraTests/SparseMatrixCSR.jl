module SparseMatrixCSRTests

    using SparseArrays
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
                CSR = sparsecsr(SparseMatrixCSR{Bi,Tv,Ti},I, J, V,maxrows,maxcols)

                @test is_entry_stored(SparseMatrixCSR{Bi,Tv,Ti},1,1)
                @test is_entry_stored(SparseMatrixCSR{Bi,Tv,Ti},1,2)
                @test is_entry_stored(SparseMatrixCSR{Bi,Tv,Ti},2,1)

                _l = 10
                _I, _J, _V = allocate_coo_vectors(SymSparseMatrixCSR{Bi,Tv,Ti},_l)
                @test length(_I) == _l
                @test length(_J) == _l
                @test length(_V) == _l
                @test eltype(_I) == Ti
                @test eltype(_J) == Ti
                @test eltype(_V) == Tv

                V2 = copy(V)
                V2 .= 0
                CSR2 = sparsecsr(SparseMatrixCSR{Bi,Tv,Ti},I, J, V2,maxrows,maxcols)
                copy_entries!(CSR2,CSR)
                @test CSR2 == CSR
                copy_entries!(CSR2,CSR2)
                @test CSR2 == CSR

                show(devnull, CSR);

                @test CSC == CSR
                @test nnz(CSC) == count(!iszero, CSC) == nnz(CSR) == count(!iszero, CSR)

                @test hasrowmajororder(CSR) == true
                @test hascolmajororder(CSR) == false
                @test getptr(CSR)           == CSR.rowptr
                @test getindices(CSR)       == colvals(CSR)

                TCSC = sparse(J, I, V, maxrows, maxcols)
                TCSR = sparsecsr(SparseMatrixCSR{Bi}, J, I, V, maxrows, maxcols)

                @test size(CSC)==size(CSR)==reverse(size(TCSC))==reverse(size(TCSC))

                @test [nzrange(CSC,col) for col in 1:size(CSC,2)] == [nzrange(TCSR,row) for row in 1:size(TCSR,1)]
                @test [nzrange(CSR,row) for row in 1:size(CSR,1)] == [nzrange(TCSC,col) for col in 1:size(TCSC,2)]

                @test nonzeros(CSC) == nonzeros(TCSR) && nonzeros(CSR) == nonzeros(TCSC)

                ICSC,JCSC,VCSC= findnz(CSC)
                ICSR,JCSR,VCSR= findnz(CSR)

                @test sort(ICSC)==sort(ICSR) && sort(JCSC)==sort(JCSR) && sort(VCSC)==sort(VCSR)

                v = rand(size(CSC)[2])
                @test CSC*v == CSR*v

                for cBi in Bi_types
                    if Bi == cBi
                        @test convert(SparseMatrixCSR{Bi}, CSR) === CSR
                    else
                        CSRC = convert(SparseMatrixCSR{cBi}, CSR)
                        @test CSRC == CSR
                        @test CSRC !== CSR
                    end
                end

                for cTi in int_types
                    for cTv in float_types
                        if (Ti,Tv) == (cTi,cTv)
                            @test convert(SparseMatrixCSR{Bi,Tv,Ti}, CSR) === CSR
                        else
                            CSRC = convert(SparseMatrixCSR{Bi,cTv,cTi}, CSR)
                            @test CSRC == CSR
                            @test CSRC !== CSR
                        end
                    end
                end

                vold = getindex(CSR,maxrows,maxcols)
                add_entry!(CSR,1,maxrows,maxcols,+)
                @test getindex(CSR,maxrows,maxcols) == vold+1

                fill_entries!(CSR,0)
                @test all(x->x==0, nonzeros(CSR))

            end

            function test_rectangular(Ti::Type,Tv::Type,nrows,ncols)
               I = Vector{Ti}()
               J = Vector{Ti}()
               V = Vector{Tv}()
               for (ik, jk, vk) in zip(rand(1:nrows, maxnz), rand(1:ncols, maxnz), rand(1:Tv(maxnz), maxnz-1))
                 push_coo!(SparseMatrixCSR{1,Tv,Ti},I,J,V,ik,jk,vk)
               end
               push_coo!(SparseMatrixCSR{1,Tv,Ti},I,J,V,nrows,ncols,maxnz)
               finalize_coo!(SparseMatrixCSR{1,Tv,Ti},I,J,V,nrows,ncols)
               CSC = sparse(I, J, V, nrows,ncols)
               CSR = sparsecsr(I, J, V, nrows,ncols)
               @test CSR == CSC
           end
           test_rectangular(Ti,Tv,maxrowsorcols,maxcols)
           test_rectangular(Ti,Tv,maxrowsorcols,maxcols)
        end
    end
end
