
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
  T = eltype(V)
  zeros(T,n)
end

"""
    allocate_in_range(::Type{V},matrix) where V

Allocate a vector of type `V` in the range of matrix `matrix`.
"""
function allocate_in_range(::Type{V},matrix) where V
  n = size(matrix,1)
  allocate_vector(V,n)
end

"""
    allocate_in_domain(::Type{V},matrix) where V

Allocate a vector of type `V` in the domain of matrix `matrix`.
"""
function allocate_in_domain(::Type{V},matrix) where V
  n = size(matrix,2)
  allocate_vector(V,n)
end

"""
    fill_entries!(a,v)

Fill the entries of array `a` with the value `v`. Returns `a`.
"""
function fill_entries!(a,v)
  fill!(a,v)
  a
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
    add_entries!(a,b,combine=+)

Perform the operation `combine` element-wise in the entries of arrays `a` and `b`
and store the result in array `a`. Returns `a`.
"""
function add_entries!(a,b,combine=+)
  @assert length(a) == length(b)
  @inbounds for i in eachindex(a)
    a[i] = combine(a[i],b[i])
  end
  a
end

"""
    add_entry!(A,v,i,j,combine=+)

Add an entry given its position and the operation to perform.
"""
function add_entry!(A::AbstractArray,v,i::Integer,j::Integer,combine=+)
  aij = A[i,j]
  A[i,j] = combine(aij,v)
end

"""
    add_entry!(A,v,i,combine=+)

Add an entry given its position and the operation to perform.
"""
function add_entry!(A::AbstractArray,v,i::Integer,combine=+)
  ai = A[i]
  A[i] = combine(ai,v)
end

"""
    scale_entries!(a,v)

Scale the entries of array `a` with the value `v`. Returns `a`.
"""
function scale_entries!(a,b)
  @inbounds for i in eachindex(a)
    a[i] = b*a[i]
  end
  a
end

# Base.mul!

"""
    muladd!(c,a,b)

Matrix multiply a*b and add to result to c. Returns c.
"""
function muladd!(c,a,b)
  _muladd!(c,a,b)
  c
end

@static if VERSION >= v"1.3"
  function _muladd!(c,a,b)
    mul!(c,a,b,1,1)
  end
else
  function _muladd!(c,a,b)
    @assert length(c) == size(a,1)
    @assert length(b) == size(a,2)
    @inbounds for j in 1:size(a,2)
      for i in 1:size(a,1)
        c[i] += a[i,j]*b[j]
      end
    end
  end
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
#
# Create the final array
# from the nz values
#
#    c = create_from_nz(b)
#
# We can also do a loop and update
# the entries of c
#
#    fill_entries!(c,0)
#    add_entry!(c,nothing,i,j)
#    add_entry!(c,v,i,j)
#

struct Loop end
struct DoNotLoop end
LoopStyle(::Type) = DoNotLoop()
LoopStyle(::T) where T = LoopStyle(T)

# For dense arrays

struct ArrayCounter{T,A}
  axes::A
  function ArrayCounter{T}(axes::A) where {T,A<:Tuple{Vararg{AbstractUnitRange}}}
    new{T,A}(axes)
  end
end

LoopStyle(::Type{<:ArrayCounter}) = DoNotLoop()

@inline add_entry!(a::ArrayCounter,args...) = nothing

nz_counter(::Type{T},axes) where T = ArrayCounter{T}(axes)

nz_allocation(a::ArrayCounter{T}) where T = fill!(similar(T,a.axes),zero(eltype(T)))

create_from_nz(a::AbstractArray) = a

# For sparse matrices

mutable struct SparseMatrixCounter{T,A}
  nnz::Int
  axes::A
  function SparseMatrixCounter{T}(axes::A) where {T,A<:NTuple{2,AbstractUnitRange}}
    nnz = 0
    new{T,A}(nnz,axes)
  end
end

LoopStyle(::Type{<:SparseMatrixCounter}) = Loop()

@inline function add_entry!(a::SparseMatrixCounter{T},::Nothing,i,j) where T
  if is_entry_stored(T,i,j)
    a.nnz = a.nnz + 1
  end
  nothing
end

@inline function add_entry!(a::SparseMatrixCounter{T},v,i,j) where T
  if is_entry_stored(T,i,j)
    a.nnz = a.nnz + 1
  end
  nothing
end

struct CooAllocation{T,A,B,C}
  counter::SparseMatrixCounter{T,A}
  I::B
  J::B
  V::C
end

LoopStyle(::Type{<:CooAllocation}) = Loop()

@inline function add_entry!(a::CooAllocation{T},::Nothing,i,j) where T
  if is_entry_stored(T,i,j)
    a.counter.nnz = a.counter.nnz + 1
    k = a.counter.nnz
    a.I[k] = i
    a.J[k] = j
  end
  nothing
end

@inline function add_entry!(a::CooAllocation{T},v,i,j) where T
  if is_entry_stored(T,i,j)
    a.counter.nnz = a.counter.nnz + 1
    k = a.counter.nnz
    a.I[k] = i
    a.J[k] = j
    a.V[k] = v
  end
  nothing
end
function nz_counter(::Type{T},axes) where T<:AbstractSparseMatrix
  SparseMatrixCounter{T}(axes)
end

function nz_allocation(a::SparseMatrixCounter{T}) where T
  counter = SparseMatrixCounter{T}(a.axes)
  I,J,V = allocate_coo_vectors(T,a.nnz)
  CooAllocation(counter,I,J,V)
end

function create_from_nz(a::CooAllocation{T}) where T
  m,n = map(length,a.counter.axes)
  finalize_coo!(T,a.I,a.J,a.V,m,n)
  sparse_from_coo(T,a.I,a.J,a.V,m,n)
end

