
"""
    rewind_ptrs!(ptrs)

Rewind the given vector of pointers.
"""
function rewind_ptrs!(ptrs::AbstractVector{<:Integer})
  @inbounds for i in (length(ptrs)-1):-1:1
    ptrs[i+1] = ptrs[i]
  end
  ptrs[1] = 1
end

"""
    length_to_ptrs!(ptrs)

Given a vector of integers, mutate it from length state to pointer state.
"""
function length_to_ptrs!(ptrs::AbstractArray{<:Integer})
  ptrs[1] = 1
  @inbounds for i in 1:(length(ptrs)-1)
    ptrs[i+1] += ptrs[i]
  end
end

function allocate_matrix end
function allocate_matrix_and_vector end

"""
    allocate_vector(::Type{V},indices) where V

Allocate a vector of type `V` indexable at the indices `indices`
"""
function allocate_vector(::Type{V},indices) where V
  n = length(indices)
  allocate_vector(V,n)
end

function allocate_vector(::Type{V},n::Integer) where V
  V(undef,n)
end

function allocate_vector(::Type{<:BlockVector{T,VV}},indices::AbstractBlockedUnitRange) where {T,VV}
  V = eltype(VV)
  mortar(map(ids -> allocate_vector(V,ids),blocks(indices)))
end

"""
    allocate_in_range(::Type{V},matrix) where V

Allocate a vector of type `V` in the range of matrix `matrix`.
"""
function allocate_in_range(::Type{V},matrix) where V
  rows = axes(matrix,1)
  allocate_vector(V,rows)
end

"""
    allocate_in_range(matrix::AbstractMatrix{T}) where T

Allocate a vector in the range of matrix `matrix`.
"""
function allocate_in_range(matrix::AbstractMatrix{T}) where T
  allocate_in_range(Vector{T},matrix)
end

function allocate_in_range(matrix::BlockMatrix{T}) where T
  V = BlockVector{T,Vector{Vector{T}}}
  allocate_in_range(V,matrix)
end

"""
    allocate_in_domain(::Type{V},matrix) where V

Allocate a vector of type `V` in the domain of matrix `matrix`.
"""
function allocate_in_domain(::Type{V},matrix) where V
  cols = axes(matrix,2)
  allocate_vector(V,cols)
end

"""
    allocate_in_domain(matrix::AbstractMatrix{T}) where T

Allocate a vector in the domain of matrix `matrix`.
"""
function allocate_in_domain(matrix::AbstractMatrix{T}) where T
  allocate_in_domain(Vector{T},matrix)
end

function allocate_in_domain(matrix::BlockMatrix{T}) where T
  V = BlockVector{T,Vector{Vector{T}}}
  allocate_in_domain(V,matrix)
end

"""
    copy_entries!(a,b)

Copy the entries of array `b` into array `a`. Returns `a`.
"""
function copy_entries!(a,b)
  if a !== b
    copyto!(a,b)
  end
  a
end

"""
    add_entry!(combine::Function,A,v,i...)
    add_entry!(A,v,i...)

Add an entry. Returns A.
"""
@inline function add_entry!(combine::Function,args...)
  @abstractmethod
end

@inline function add_entry!(args...)
  add_entry!(+,args...)
end

@inline function add_entry!(combine::Function,A::AbstractMatrix,v,i,j)
  aij = A[i,j]
  A[i,j] = combine(aij,v)
  A
end

@inline function add_entry!(combine::Function,A::AbstractVector,v,i)
  ai = A[i]
  A[i] = combine(ai,v)
  A
end

@inline function add_entry!(combine::Function,A::AbstractMatrix,v::Nothing,i,j)
  A
end

@inline function add_entry!(combine::Function,A::AbstractVector,v::Nothing,i)
  A
end

# Warning: the usage of @inline and @noinline seems to have dramatic performance
# implications. Do not change it.

@noinline function add_entries!(A,vs,is,js)
  add_entries!(+,A,vs,is,js)
end

@noinline function add_entries!(A,vs,is)
  add_entries!(+,A,vs,is)
end

"""
    add_entries!(combine::Function,A,vs,is,js)

Add several entries only for positive input indices. Returns A.
"""
@inline function add_entries!(combine::Function,A,vs,is,js)
  _add_entries!(combine,A,vs,is,js)
end

@inline function _add_entries!(combine::Function,A,vs::Nothing,is,js)
  for (lj,j) in enumerate(js)
    if j>0
      for (li,i) in enumerate(is)
        if i>0
          add_entry!(combine,A,nothing,i,j)
        end
      end
    end
  end
  A
end

@inline function _add_entries!(combine::Function,A,vs,is,js)
  for (lj,j) in enumerate(js)
    if j>0
      for (li,i) in enumerate(is)
        if i>0
          vij = vs[li,lj]
          add_entry!(combine,A,vij,i,j)
        end
      end
    end
  end
  A
end

@inline function add_entries!(combine::Function,A,vs,is)
  _add_entries!(combine,A,vs,is)
end

@inline function _add_entries!(combine::Function,A,vs::Nothing,is)
  for (li, i) in enumerate(is)
    if i>0
      add_entry!(A,nothing,i)
    end
  end
  A
end

@inline function _add_entries!(combine::Function,A,vs,is)
  for (li, i) in enumerate(is)
    if i>0
      vi = vs[li]
      add_entry!(A,vi,i)
    end
  end
  A
end

@inline function add_entries!(combine::Function,A::AbstractMatrix,vs::Nothing,is,js)
  A
end

@inline function add_entries!(combine::Function,A::AbstractVector,vs::Nothing,is)
  A
end


"""
    muladd!(c,a,b)

Matrix multiply a*b and add to result to c. Returns c.
"""
muladd!(b,A,x) = mul!(b,A,x,one(eltype(b)),one(eltype(b)))

"""
    axpy_entries!(α::Number, A::T, B::T) where {T<: AbstractMatrix} -> T

Efficient implementation of axpy! for sparse matrices.
"""
function axpy_entries!(α::Number, A::T, B::T) where {T<:AbstractMatrix}
  iszero(α) && return B

  axpy!(α, A, B)
  B
end

# For sparse matrices, it is surprisingly quicker to call `@. B += α * A` than
# `axpy!(α, A, B)`.` Calling axpy! on the nonzero values of A and B is the most
# efficient approach but this is only possible when A and B have the same
# sparsity pattern. The checks add some non-negligible overhead so we make them
# optional by adding a keyword.
const cannot_axpy_entries_msg = """
It is only possible to efficiently add two sparse matrices that have the same
sparsity pattern.
"""

function axpy_entries!(
  α::Number, A::T, B::T;
  check::Bool=true
) where {T<:SparseMatrixCSC}
  iszero(α) && return B

  if check
    msg = cannot_axpy_entries_msg
    @check rowvals(A) == rowvals(B) msg
    @check all(nzrange(A, j) == nzrange(B, j) for j in axes(A, 2)) msg
  end

  axpy!(α, nonzeros(A), nonzeros(B))
  B
end

function axpy_entries!(
  α::Number, A::T, B::T;
  check::Bool=true
) where {T<:Union{SparseMatrixCSR,SymSparseMatrixCSR}}
  iszero(α) && return B

  if check
    msg = cannot_axpy_entries_msg
    @check colvals(A) == colvals(B) msg
    @check all(nzrange(A, j) == nzrange(B, j) for j in axes(A, 1)) msg
  end

  axpy!(α, nonzeros(A), nonzeros(B))
  B
end

function axpy_entries!(α::Number, A::T, B::T) where {T<:AbstractBlockMatrix}
  map(blocks(A), blocks(B)) do a, b
    axpy_entries!(α, a, b)
  end
  B
end

#
# Some API associated with assembly routines
#
# Generate a counter to count the nz values
# for an array type A
#
#    a = nz_counter(A,(rows,cols))
#
# Do a loop to count the number of nz values
# Do the loop only when needed by using the LoopStyle
# For instance, when assembling a dense vector the loop
# is not needed
# A nz value is counted by calling add_entry!
#
#    if LoopStyle(a) == Loop()
#      add_entry!(a,nothing,i,j)
#      add_entry!(a,v,i,j)
#      add_entries!(a,nothing,is,js)
#      add_entries!(a,vs,is,js)
#    end
#
# Now we can allocate the nz values
# This can be the vectors in coo format or
# it can already be the final array for dense vectors
#
#    b = nz_allocation(a)
#
# Do a loop to set entries
# We can use nothing if we want to add an entry to the sparsity
# pattern but we don't want to set a value to it yet
#
#    add_entry!(b,nothing,i,j)
#    add_entry!(b,v,i,j)
#    add_entries!(b,nothing,is,js)
#    add_entries!(b,vs,is,js)
#
# Create the final array
# from the nz values
#
#    c = create_from_nz(b)
#
# We can also do a loop and update
# the entries of c
#
#    fill!(c,0) or LinearAlgebra.fillstored!(c,0)
#    add_entry!(c,v,i,j)
#    add_entries!(c,vs,is,js)
#

"""
    struct Loop end
"""
struct Loop end
"""
    struct DoNotLoop end
"""
struct DoNotLoop end
"""
    LoopStyle(::Type)
    LoopStyle(::T) where T = LoopStyle(T)

Trait to tell if looping on nonzeros entries is necessary to count the values
with a counter (see [`nz_counter`](@ref)) of an array of type `T`.

Returns [`Loop()`](@ref Loop) or [`DoNotLoop`](@ref).
"""
LoopStyle(::Type) = DoNotLoop()
LoopStyle(::T) where T = LoopStyle(T)

# By default process, matrix and vector separately
# but, in some situations, create_from_nz of the vector
# can reuse data from the one computed in
# create_from_nz for the matrix (e.g., GridapDistributed)
function create_from_nz(a,b)
  A = create_from_nz(a)
  B = create_from_nz(b)
  A,B
end
# See comment above for create_from_nz. The same applies here
# for nz_allocation.
function nz_allocation(a,b)
  nz_allocation(a),nz_allocation(b)
end

# For dense arrays

"""
    struct ArrayBuilder{T}

`T` is the type of array to be built.
"""
struct ArrayBuilder{T}
  array_type::Type{T}
end

ArrayBuilder(a::ArrayBuilder) = a

"""
    get_array_type(::ArrayBuilder{T}) = T
"""
get_array_type(::ArrayBuilder{T}) where T = T

"""
    ArrayCounter{T}(axes::A)

where `T` is the target matrix/array type.
"""
struct ArrayCounter{T,A}
  axes::A

  function ArrayCounter{T}(axes::A) where {T,A<:Tuple{Vararg{AbstractUnitRange}}}
    new{T,A}(axes)
  end
end

LoopStyle(::Type{<:ArrayCounter}) = DoNotLoop()

@inline add_entry!(c::Function,a::ArrayCounter,v,i,j) = a
@inline add_entry!(c::Function,a::ArrayCounter,v,i) = a
@inline add_entries!(c::Function,a::ArrayCounter,v,i,j) = a
@inline add_entries!(c::Function,a::ArrayCounter,v,i) = a

#nz_counter(::Type{T},axes) where T = ArrayCounter{T}(axes)

"""
    nz_counter(::ArrayBuilder{T},axes)

Generate a counter ([`ArrayCounter`](@ref)) to count the nonzero values for an array type A.
"""
nz_counter(::ArrayBuilder{T},axes) where T = ArrayCounter{T}(axes)

"""
    nz_allocation(a::ArrayCounter{T})

Allocates a vector that will serve as structural nonzero values internal storage
for a matrix holding the values counted in `a`. See also
[`create_from_nz`](@ref), `SparseArrays.nonzeros`.
"""
nz_allocation(a::ArrayCounter{T}) where T = fill!(similar(T,map(length,a.axes)),zero(eltype(T)))

"""
    create_from_nz(nz_alloc::AbstractArray)

Creates a matrix from its nonzero values, allocated using [`nz_allocation`](@ref)
and filled using [`add_entry!`](@ref) and/or [`add_entries!`](@ref).
"""
create_from_nz(a::AbstractArray) = a

# For sparse matrices

"""
    struct MinCPU end

Represent sparse matrix assembly strategy minimizing CPU cost.
"""
struct MinCPU end

"""
    struct MinMemory{T}

Represent sparse matrix assembly strategy minimizing memory cost.
"""
struct MinMemory{T}
  maxnnz::T
end

MinMemory() = MinMemory(nothing)

"""
    struct SparseMatrixBuilder{T,A}
"""
struct SparseMatrixBuilder{T,A}
  matrix_type::Type{T}
  approach::A
end

"""
    SparseMatrixBuilder(::Type{T}[, ::S])

Create a builder for sparse matrix of type `T`, using assembly stategy `S()`,
either [`MinCPU()`](@ref MinCPU) or [`MinMemory()`](@ref MinMemory).
The default is `MinMemory()`.
"""
SparseMatrixBuilder(::Type{T}) where T = SparseMatrixBuilder(T,MinMemory())
SparseMatrixBuilder(a::SparseMatrixBuilder) = a

get_array_type(::SparseMatrixBuilder{T}) where T = T

mutable struct CounterCOO{T,A}
  nnz::Int
  axes::A
  function CounterCOO{T}(axes::A) where {T,A<:NTuple{2,AbstractUnitRange}}
    nnz = 0
    new{T,A}(nnz,axes)
  end
end

LoopStyle(::Type{<:CounterCOO}) = Loop()

@inline function add_entry!(::Function,a::CounterCOO{T},v,i,j) where T
  if is_entry_stored(T,i,j)
    a.nnz = a.nnz + 1
  end
  a
end

struct AllocationCOO{T,A,B,C}
  counter::CounterCOO{T,A}
  I::B
  J::B
  V::C
end

LoopStyle(::Type{<:AllocationCOO}) = Loop()

@inline function add_entry!(::typeof(+),a::AllocationCOO{T},::Nothing,i,j) where T
  if is_entry_stored(T,i,j)
    a.counter.nnz = a.counter.nnz + 1
    k = a.counter.nnz
    a.I[k] = i
    a.J[k] = j
  end
  nothing
end

@inline function add_entry!(::typeof(+),a::AllocationCOO{T},v,i,j) where T
  if is_entry_stored(T,i,j)
    a.counter.nnz = a.counter.nnz + 1
    k = a.counter.nnz
    a.I[k] = i
    a.J[k] = j
    a.V[k] = v
  end
  nothing
end

#function nz_counter(::Type{T},axes) where T<:AbstractSparseMatrix
#  CounterCOO{T}(axes)
#end

function nz_counter(::SparseMatrixBuilder{T},axes) where T<:AbstractSparseMatrix
  CounterCOO{T}(axes)
end

function nz_allocation(a::CounterCOO{T}) where T
  counter = CounterCOO{T}(a.axes)
  I,J,V = allocate_coo_vectors(T,a.nnz)
  AllocationCOO(counter,I,J,V)
end

function create_from_nz(a::AllocationCOO{T}) where T
  m,n = map(length,a.counter.axes)
  finalize_coo!(T,a.I,a.J,a.V,m,n)
  sparse_from_coo(T,a.I,a.J,a.V,m,n)
end

# The following methods can be implemented for sparse matrices
# instead of nz_counter, nz_allocation, and create_from_nz

"""
    sparse_from_coo(::Type,I,J,V,m,n)
"""
function sparse_from_coo(::Type,I,J,V,m,n)
  @abstractmethod
end

"""
    is_entry_stored(::Type,i,j) -> Bool

Tells if the entry with coordinates `[i,j]` will be stored in the coo vectors.
"""
function is_entry_stored(::Type,i,j)
  @abstractmethod
end

"""
    finalize_coo!(::Type,I,J,V,m,n)

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(::Type,I,J,V,m,n)
  @abstractmethod
end

"""
    nz_index(A::AbstractSparseMatrix,i,j)

Index of `A`[`i`,`j`] in the structural nonzero values internal storage of `A`.
See also `SparseArrays.nonzeros`.
"""
function nz_index(A::AbstractSparseMatrix,i,j)
  @abstractmethod
end

"""
    push_coo!(::Type, I,J,V,i,j,v)
Inserts entries in COO vectors for further building a sparse matrix of type T.
"""
function push_coo!(::Type,I,J,V,i,j,v)
  @abstractmethod
end

function add_entry!(combine::Function,A::AbstractSparseMatrix,v::Number,i,j)
  k = nz_index(A,i,j)
  nz = nonzeros(A)
  Aij = nz[k]
  nz[k] = combine(Aij,v)
  A
end

function copy_entries!(a::T,b::T) where T<:AbstractSparseMatrix
  na = nonzeros(a)
  nb = nonzeros(b)
  if na !== nb
    copyto!(na,nb)
  end
end

"""
    allocate_coo_vectors(::Type{<:AbstractSparseMatrix{Tv,Ti}}, n::Integer)
"""
function allocate_coo_vectors(
   ::Type{<:AbstractSparseMatrix{Tv,Ti}},n::Integer) where {Tv,Ti}
  (zeros(Ti,n), zeros(Ti,n), zeros(Tv,n))
end
