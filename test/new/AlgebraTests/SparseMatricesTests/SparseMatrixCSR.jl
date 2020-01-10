@testset "SparseMatrixCSR" begin
    maxnz=5
    maxrows=5
    maxcols=5
    int_types=(Int32,Int64)
    float_types=(Float32,Float64)
    Bi_types=(0,1)

    for Ti in int_types
        for Tv in float_types
            for Bi in Bi_types
                I = Vector{Ti}()
                J = Vector{Ti}()
                V = Vector{Tv}()
                for (ik, jk, vk) in zip(rand(1:maxrows, maxnz), rand(1:maxcols, maxnz), rand(1:Tv(maxnz), maxnz))
                    push_coo!(SparseMatrixCSR{Bi,Tv,Ti},I,J,V,ik,jk,vk)
                end
                finalize_coo!(SparseMatrixCSR{Bi,Tv,Ti},I,J,V,maxrows,maxcols)
                CSC = sparse(I, J, V, maxrows,maxcols)
                CSR = sparsecsr(SparseMatrixCSR{Bi,Tv,Ti},I, J, V,maxrows,maxcols)

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

                @test sort(ICSC)==sort(JCSR) && sort(JCSC)==sort(ICSR) && sort(VCSC)==sort(VCSR)

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

            end
        end
    end    
end

