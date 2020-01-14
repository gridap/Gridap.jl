@testset "SymSparseMatrixCSR" begin
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
                for (ik, jk, vk) in zip(rand(1:maxrows, maxnz), rand(1:maxcols, maxnz), rand(1:Tv(maxnz), maxnz-1) )
                    push_coo!(SymSparseMatrixCSR{Bi,Tv,Ti},I,J,V,ik,jk,vk)
                end
                push_coo!(SymSparseMatrixCSR{Bi,Tv,Ti},I,J,V,maxrows,maxcols,maxnz)
                finalize_coo!(SymSparseMatrixCSR{Bi,Tv,Ti},I,J,V,maxrows, maxcols)
                SYMCSC = Symmetric(sparse(I, J, V, maxrows, maxcols),:U)
                SYMCSR = symsparsecsr(SymSparseMatrixCSR{Bi,Tv,Ti},I, J, V, maxrows, maxcols)

                @test size(SYMCSC)==size(SYMCSR)
                @test SYMCSC == SYMCSR

                @test convert(SymSparseMatrixCSR{Bi}, SYMCSR) === SYMCSR

                @test hasrowmajororder(SYMCSR) == true
                @test hascolmajororder(SYMCSR) == false
                @test getptr(SYMCSR)           == SYMCSR.uppertrian.rowptr
                @test getindices(SYMCSR)       == colvals(SYMCSR)

                @test nnz(SYMCSC.data) == nnz(SYMCSR.uppertrian) <= nnz(SYMCSR) 
                @test count(!iszero, SYMCSC.data) == count(!iszero, SYMCSR.uppertrian)

                ICSC,JCSC,VCSC= findnz(SYMCSC.data)
                ICSR,JCSR,VCSR= findnz(SYMCSR)

                @test sort(ICSC)==sort(JCSR) && sort(JCSC)==sort(ICSR) && sort(VCSC)==sort(VCSR)

                v = rand(size(SYMCSC)[2])
                @test SYMCSC*v == SYMCSR*v

                for cBi in Bi_types
                    if Bi == cBi
                        @test convert(SymSparseMatrixCSR{Bi,Tv,Ti}, SYMCSR) === SYMCSR
                    else
                        SYMCSRC = convert(SymSparseMatrixCSR{cBi,Tv,Ti}, SYMCSR)
                        @test SYMCSRC == SYMCSR
                        @test SYMCSRC !== SYMCSR
                    end
                end

                for cTi in int_types
                    for cTv in float_types
                        if (Ti,Tv) == (cTi,cTv)
                            @test convert(SymSparseMatrixCSR{Bi,Tv,Ti}, SYMCSR) === SYMCSR
                        else
                            SYMCSRC = convert(SymSparseMatrixCSR{Bi,cTv,cTi}, SYMCSR)
                            @test SYMCSRC == SYMCSR
                            @test SYMCSRC !== SYMCSR
                        end
                    end
                end

                vold = getindex(SYMCSR,maxrows,maxcols)
                add_entry!(SYMCSR,1,maxrows,maxcols,+)
                @test getindex(SYMCSR,maxrows,maxcols) == vold+1

                fill_entries!(SYMCSR,0)
                @test all(x->x==0, nonzeros(SYMCSR))

            end
        end
    end
    
end

