
"""
    struct SparseMatrixCSR{Bi,Tv,Ti<:Integer} <: AbstractSparseMatrix{Tv,Ti}

Matrix type for storing Bi-based sparse matrices 
in the Compressed Sparse Row format. The standard 
way of constructing SparseMatrixCSR is through the 
[`sparsecsr`](@ref) function.
"""
struct SparseMatrixCSR{Bi,Tv,Ti} <: AbstractSparseMatrix{Tv,Ti}
    offset :: Int             # Bi-index offset (Bi-1)
    m      :: Int             # Number of columns
    n      :: Int             # Number of rows
    rowptr :: Vector{Ti}      # Row i is in (rowptr[i]+offset):(rowptr[i+1]-Bi)
    colval :: Vector{Ti}      # Col indices (with base Bi) of stored values
    nzval  :: Vector{Tv}      # Stored values, typically nonzeros

    function SparseMatrixCSR{Bi}(m::Integer, n::Integer, rowptr::Vector{Ti}, colval::Vector{Ti},
                                    nzval::Vector{Tv}) where {Bi,Tv,Ti<:Integer}
        @noinline throwsz(str, lbl, k) =
            throw(ArgumentError("number of $str ($lbl) must be ≥ 0, got $k"))
        m < 0 && throwsz("rows", 'm', m)
        n < 0 && throwsz("columns", 'n', n)
        rowptr[m+1] == length(nzval) && throwsz("nnz", "rowptr[m+1]", length(nzval))
        offset = Bi-1
        rowptr .+= offset
        colval .+= offset
        new{Bi,Tv,Ti}(Int(offset), Int(m), Int(n), rowptr, colval, nzval)
    end
end


SparseMatrixCSR(transpose::SparseMatrixCSC) where {Tv,Ti} = 
        SparseMatrixCSR{1}(transpose.n, transpose.m, transpose.colptr, transpose.rowval, transpose.nzval)
SparseMatrixCSR{Bi}(transpose::SparseMatrixCSC{Tv,Ti}) where {Bi,Tv,Ti} = 
        SparseMatrixCSR{Bi}(transpose.n, transpose.m, transpose.colptr, transpose.rowval, transpose.nzval)

size(S::SparseMatrixCSR) = (S.m, S.n)
show(io::IO, S::SparseMatrixCSR) = show(convert(IOContext, io), S::SparseMatrixCSR)

hasrowmajororder(::Type{<:SparseMatrixCSR}) = true
hasrowmajororder(a::SparseMatrixCSR) = hasrowmajororder(SparseMatrixCSR)
hascolmajororder(::Type{<:SparseMatrixCSR}) = false
hascolmajororder(a::SparseMatrixCSR) = hascolmajororder(SparseMatrixCSR)

getptr(S::SparseMatrixCSR) = S.rowptr
getindices(S::SparseMatrixCSR) = colvals(S)
nnz(S::SparseMatrixCSR{Bi}) where {Bi} = length(nonzeros(S))
count(pred, S::SparseMatrixCSR) = count(pred, nonzeros(S))
nonzeros(S::SparseMatrixCSR) = S.nzval
nzrange(S::SparseMatrixCSR{Bi}, row::Integer) where {Bi} = S.rowptr[row]-S.offset:S.rowptr[row+1]-Bi
rowvals(S::SparseMatrixCSR) = error("CSR sparse matrix does not contain raw row values")
colvals(S::SparseMatrixCSR) = S.colval

sparsecsr(I,J,args...) = 
        SparseMatrixCSR(sparse(J,I,args...))
sparsecsr(::Type{<:SparseMatrixCSR},I,J,args...) = 
        sparsecsr(I,J,args...)
sparsecsr(::Type{<:SparseMatrixCSR{Bi}},I,J,args...) where {Bi} = 
        SparseMatrixCSR{Bi}(sparse(J,I,args...))
sparsecsr(::Type{<:SparseMatrixCSR{Bi}},I,J,V,m,n,args...) where {Bi} = 
        SparseMatrixCSR{Bi}(sparse(J,I,V,n,m,args...))
sparsecsr(::Type{<:SparseMatrixCSR{Bi,Tv,Ti}},I::Vector{Ti},J::Vector{Ti},args...) where {Bi,Tv,Ti} = 
        SparseMatrixCSR{Bi}(sparse(J,I,args...))
sparsecsr(::Type{<:SparseMatrixCSR{Bi,Tv,Ti}},I::Vector{Ti},J::Vector{Ti},V::Vector{Tv},m,n,args...) where {Bi,Tv,Ti} = 
        SparseMatrixCSR{Bi}(sparse(J,I,V,n,m,args...))

*(A::SparseMatrixCSR, v::Vector) = (y = similar(v,A.n);mul!(y,A,v))

convert(::Type{<:SparseMatrixCSR}, x::AbstractSparseMatrix)  = convert(SparseMatrixCSR{1}, x)

function show(io::IO, ::MIME"text/plain", S::SparseMatrixCSR)
    xnnz = nnz(S)
    print(io, S.m, "×", S.n, " ", typeof(S), " with ", xnnz, " stored ",
              xnnz == 1 ? "entry" : "entries")
    if xnnz != 0
        print(io, ":")
        show(IOContext(io, :typeinfo => eltype(S)), S)
    end
end


function show(io::IOContext, S::SparseMatrixCSR{Bi}) where{Bi}
    nnz(S) == 0 && return show(io, MIME("text/plain"), S)

    ioc = IOContext(io, :compact => true)
    function _format_line(r, col, padr, padc)
        print(ioc, "\n  [", rpad(col, padr), ", ", lpad(S.colval[r], padc), "]  =  ")
        if isassigned(S.nzval, Int(r))
            show(ioc, S.nzval[r])
        else
            print(ioc, Base.undef_ref_str)
        end
    end

    function _get_cols(from, to)
        idx = eltype(S.rowptr)[]
        c = searchsortedlast(S.rowptr, from)
        for i = from:to
            while i == S.rowptr[c+1]
                c +=1
            end
            push!(idx, c)
        end
        idx
    end

    rows = displaysize(io)[1] - 4 # -4 from [Prompt, header, newline after elements, new prompt]
    if !get(io, :limit, false) || rows >= nnz(S) # Will the whole matrix fit when printed?
        cols = _get_cols(Bi, nnz(S)+S.offset)
        padr, padc = ndigits.((maximum(S.colval[1:nnz(S)]), cols[end]))
        _format_line.(1:nnz(S), cols.+S.offset, padr, padc)
    else
        if rows <= 2
            print(io, "\n  \u22ee")
            return
        end
        s1, e1 = 1, div(rows - 1, 2) # -1 accounts for \vdots
        s2, e2 = nnz(S) - (rows - 1 - e1) + 1, nnz(S)
        cols1, cols2 = _get_cols(s1+S.offset, e1+S.offset), _get_cols(s2+S.offset, e2+S.offset)
        padr = ndigits(max(maximum(S.colval[s1:e1]), maximum(S.colval[s2:e2])))
        padc = ndigits(cols2[end])
        _format_line.(s1:e1, cols1, padr, padc)
        print(io, "\n  \u22ee")
        _format_line.(s2:e2, cols2, padr, padc)
    end
    return
end


function getindex(S::SparseMatrixCSR{Bi,Tv,Ti}, i0::Integer, i1::Integer) where {Bi,Tv,Ti}
    if !(1 <= i0 <= S.m && 1 <= i1 <= S.n); throw(BoundsError()); end
    r1 = Int(S.rowptr[i0]-S.offset)
    r2 = Int(S.rowptr[i0+1]-Bi)
    (r1 > r2) && return zero(Tv)
    i1_= i1+S.offset
    r1 = searchsortedfirst(S.colval, i1_, r1, r2, Base.Order.Forward)
    ((r1 > r2) || (S.colval[r1] != i1_)) ? zero(Tv) : S.nzval[r1]
end


function findnz(S::SparseMatrixCSR{Bi,Tv,Ti}) where {Bi,Tv,Ti}
    numnz = nnz(S)
    I = Vector{Ti}(undef, numnz)
    J = Vector{Ti}(undef, numnz)
    V = Vector{Tv}(undef, numnz)

    count = 1
    @inbounds for row in 1:S.m
        @inbounds for k in nzrange(S,row)
            I[count] = S.colval[k]-S.offset
            J[count] = row
            V[count] = S.nzval[k]
            count += 1
        end
    end

    return (I, J, V)
end

function sparse_from_coo(M::Type{<:SparseMatrixCSR}, args...)
  sparsecsr(M, args...)
end

function push_coo!(::Type{<:SparseMatrixCSR},
    I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number) where {Bi}
    (push!(I, ik), push!(J, jk), push!(V, vk))
end

function finalize_coo!(::Type{<:SparseMatrixCSR},
    I::Vector,J::Vector,V::Vector,m::Integer,n::Integer) 
end

function add_entry!(A::SparseMatrixCSR{Bi,Tv,Ti},v::Number,i::Integer,j::Integer,combine::Function=+) where {Bi,Tv,Ti<:Integer}
    cv = convert(Tv, v)
    ci = convert(Ti, i)
    cj = convert(Ti, j+A.offset)
    m,n = size(A)
    if !((1 <= i <= m) & (1 <= j <= n))
        throw(BoundsError(A, (i,j)))
    end
    rows = getptr(A)
    indices = getindices(A)
    vals = nonzeros(A)
    rowifirstk = Int(rows[ci]-A.offset)
    rowilastk = Int(rows[ci+1]-Bi)
    searchk = searchsortedfirst(indices, cj, rowifirstk, rowilastk, Base.Order.Forward)
    if searchk <= rowilastk-A.offset && indices[searchk] == cj
        # Column j contains entry A[i,j]. Update
        vold = vals[searchk]
        vals[searchk] = combine(vold, cv)
    else
        throw(ArgumentError("Entry ($i,$j) does not exists in the sparsity pattern"))
    end
    return A
end


function mul!(y::AbstractVector,A::SparseMatrixCSR,v::AbstractVector{T}) where {T}
    A.n == size(v, 1) || throw(DimensionMismatch())
    A.m == size(y, 1) || throw(DimensionMismatch())

    y .= zero(T)
    for row = 1:size(y, 1)
        @inbounds for nz in nzrange(A,row)
            col = A.colval[nz]-A.offset
            y[row] += A.nzval[nz]*v[col]
        end
    end
    return y
end


function convert(::Type{<:SparseMatrixCSR{Bi}}, x::SparseMatrixCSR{xBi}) where {Bi,xBi}
    if Bi == xBi
        return x
    else
        return SparseMatrixCSR{Bi}( x.m, 
                                    x.n, 
                                    copy(getptr(x)).-x.offset, 
                                    copy(getindices(x)).-x.offset, 
                                    copy(nonzeros(x)))
    end
end


function convert(::Type{<:SparseMatrixCSR{Bi,Tv,Ti}}, x::SparseMatrixCSR{xBi,xTv,xTi}) where {Bi,Tv,Ti,xBi,xTv,xTi}
    if (Bi,Tv,Ti) == (xBi,xTv,xTi)
        return x
    else
        return SparseMatrixCSR{Bi}( x.m, 
                                    x.n, 
                                    convert(Vector{Ti}, copy(getptr(x)).-x.offset), 
                                    convert(Vector{Ti}, copy(getindices(x)).-x.offset), 
                                    convert(Vector{Tv}, copy(nonzeros(x))))
    end
end


function convert(::Type{<:SparseMatrixCSR{Bi}}, x::SparseMatrixCSC) where {Bi}
    A = sparse(transpose(x))
    (m, n) = size(A)
    return SparseMatrixCSR{Bi}(m, n, getptr(A), getindices(A), nonzeros(A))
end


function convert(::Type{<:SparseMatrixCSR{Bi,Tv,Ti}}, x::SparseMatrixCSC) where {Bi,Tv,Ti}
    A = sparse(transpose(x))
    return SparseMatrixCSR{Bi}( A.m, 
                                A.n, 
                                convert(Vector{Ti}, getptr(A)), 
                                convert(Vector{Ti}, getindices(A)), 
                                convert(Vector{Tv}, nonzeros(A)))
end


function convert(::Type{SparseMatrixCSC{Tv,Ti}}, x::SparseMatrixCSR{xBi,xTv,xTi}) where {Tv,Ti,xBi,xTv,xTi}
    A = sparse(transpose(x))
    return SparseMatrixCSR{Bi}( A.m, 
                                A.n, 
                                convert(Vector{Ti}, getptr(A)), 
                                convert(Vector{Ti}, getindices(A)), 
                                convert(Vector{Tv}, nonzeros(A)))
end

