
"""
    struct SymSparseMatrixCSR{T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}

Matrix type for storing symmetric sparse matrices in the
Compressed Sparse Row format. The standard way
of constructing SparseMatrixCSR is through the 
[`symsparsecsr`](@ref) function.
"""
struct SymSparseMatrixCSR{Bi,T,Ti<:Integer} <: AbstractSparseMatrix{T,Ti}
    uppertrian :: SparseMatrixCSR{Bi,T,Ti}
end

# SparseMatrix interface implementation
function push_coo!(::Type{<:SymSparseMatrixCSR},
        I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number)
    (ik>jk) && return
    (push!(I, ik), push!(J, jk), push!(V, vk))
end

function finalize_coo!(T::Type{<:SymSparseMatrixCSR},
        I::Vector,J::Vector,V::Vector, m::Integer, n::Integer) 
    m == n || throw(DimensionMismatch("matrix is not square: dimensions are ($m, $n)"))
    touched = zeros(Bool,m)
    for k in 1:length(I)
        Ik = I[k]
        Jk = J[k]
        if Ik == Jk
            touched[Ik] = true
        end
    end
    for k in 1:m
        if ! touched[k]
            push_coo!(T,I,J,V,k,k,zero(eltype(V)))
        end
    end
end

add_entry!(A::SymSparseMatrixCSR{Bi,Tv,Ti},v::Number,i::Integer,j::Integer,combine::Function=+) where {Bi,Tv,Ti<:Integer} =
        return i>j ? A : add_entry!(A.uppertrian,v,i,j,combine)

# CompressedSparseMatrix interface implementation

hasrowmajororder(::Type{<:SymSparseMatrixCSR}) = true
hasrowmajororder(a::SymSparseMatrixCSR) = hasrowmajororder(SymSparseMatrixCSR)
hascolmajororder(::Type{<:SymSparseMatrixCSR}) = false
hascolmajororder(a::SymSparseMatrixCSR) = hascolmajororder(SymSparseMatrixCSR)
getptr(S::SymSparseMatrixCSR) = getptr(S.uppertrian)
getindices(S::SymSparseMatrixCSR) = colvals(S)
getindex(A::SymSparseMatrixCSR, x::Integer, y::Integer) = getindex(A.uppertrian,min(x,y),max(x,y))
nnz(S::SymSparseMatrixCSR) = nnz(S.uppertrian)
count(pred, S::SymSparseMatrixCSR) = count(pred, S.uppertrian)
nonzeros(S::SymSparseMatrixCSR) = nonzeros(S.uppertrian)
nzrange(S::SymSparseMatrixCSR, row::Integer) = nzrange(S.uppertrian, row)
findnz(S::SymSparseMatrixCSR) = findnz(S.uppertrian)
rowvals(S::SymSparseMatrixCSR) = rowvals(S.uppertrian)
colvals(S::SymSparseMatrixCSR) = colvals(S.uppertrian)

show(io::IO, S::SymSparseMatrixCSR) = show(io, S.uppertrian)

size(A::SymSparseMatrixCSR) = size(A.uppertrian)

*(A::SymSparseMatrixCSR, v::Vector) = (y = similar(v,A.uppertrian.n);mul!(y,A,v))

function symsparsecsr(T::Type{<:SymSparseMatrixCSR{Bi,Tv,Ti}},I::Vector{Ti},J::Vector{Ti},V::Vector{Tv},args...)  where {Bi,Tv,Ti}
    m = length(args)>0 ? args[1] : isempty(I) ? 0 : Int(maximum(I))
    n = length(args)>1 ? args[2] : m
    c = length(args)>2 ? args[3] : +
    m == n || throw(DimensionMismatch("matrix is not square: dimensions are ($m, $n)"))
    SymSparseMatrixCSR(sparsecsr(SparseMatrixCSR{Bi},I,J,V,m,n,c))
end

symsparsecsr(I,J,V,args...) =
           SymSparseMatrixCSR(sparsecsr(I,J,V,args...))
symsparsecsr(T::Type{<:SymSparseMatrixCSR},I::Vector{Ti},J::Vector{Ti},V::Vector{Tv},args...) where {Tv,Ti} =
           symsparsecsr(SparseMatrixCSR{1,Tv,Ti},I,J,V,args...)
symsparsecsr(T::Type{<:SymSparseMatrixCSR{Bi}},I::Vector{Ti},J::Vector{Ti},V::Vector{Tv},args...) where {Bi,Tv,Ti} =
           symsparsecsr(SparseMatrixCSR{Bi,Tv,Ti},I,J,V,args...)


function show(io::IO, ::MIME"text/plain", S::SymSparseMatrixCSR)
    xnnz = nnz(S)
    print(io, S.uppertrian.m, "×", S.uppertrian.n, " ", 
              typeof(S), " with ", xnnz, " stored ",
              xnnz == 1 ? "entry" : "entries")
    if xnnz != 0
        print(io, ":")
        show(IOContext(io, :typeinfo => eltype(S)), S)
    end
end


function mul!(y::AbstractVector,A::SymSparseMatrixCSR,v::AbstractVector{T}) where {T}
    A.uppertrian.n == size(v, 1) || throw(DimensionMismatch())
    A.uppertrian.m == size(y, 1) || throw(DimensionMismatch())

    y .= zero(T)
    for row = 1:size(y, 1)
        @inbounds for nz in nzrange(A,row)
            col = A.uppertrian.colval[nz]-A.uppertrian.offset
            y[row] += A.uppertrian.nzval[nz]*v[col]
            row != col && (y[col] += A.uppertrian.nzval[nz]*v[row])
        end
    end
    return y
end


function convert(::Type{<:SymSparseMatrixCSR{Bi,Tv,Ti}}, x::SymSparseMatrixCSR{xBi,xTv,xTi}) where {Bi,Tv,Ti,xBi,xTv,xTi}
    if (Bi,Tv,Ti) == (xBi,xTv,xTi)
        return x
    else
        return SymSparseMatrixCSR(convert(SparseMatrixCSR{Bi,Tv,Ti}, x.uppertrian))
    end
end

