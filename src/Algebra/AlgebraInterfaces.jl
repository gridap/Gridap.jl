
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
For sparse matrices it only fills the non-zero entries.
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

"""
    add_entries!(combine::Function,A,vs,is...)

Add several entries only for positive input indices. Returns A.
"""
@inline function add_entries!(combine::Function,args...)
  @abstractmethod
end

@inline function add_entries!(args...)
  add_entries!(+,args...)
end

@inline function add_entries!(combine::Function,A,vs,is,js)
  @inline _vij(vs,i,j) = vs[i,j]
  @inline _vij(vs::Nothing,i,j) = vs
  for (lj,j) in enumerate(js)
    if j>0
      for (li,i) in enumerate(is)
        if i>0
          vij = _vij(vs,li,lj)
          add_entry!(combine,A,vij,i,j)
        end
      end
    end
  end
  A
end


@inline function add_entries!(combine::Function,A,vs,is)
  @inline _vi(vs,i) = vs[i]
  @inline _vi(vs::Nothing,i) = vs
  for (li, i) in enumerate(is)
    if i>0
      vi = _vi(vs,li)
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
#    fill_entries!(c,0)
#    add_entry!(c,v,i,j)
#    add_entries!(c,vs,is,js)
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

@inline add_entry!(c::Function,a::ArrayCounter,args...) = a
@inline add_entries!(c::Function,a::ArrayCounter,args...) = a

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

@inline function add_entry!(::Function,a::SparseMatrixCounter{T},v,i,j) where T
  if is_entry_stored(T,i,j)
    a.nnz = a.nnz + 1
  end
  a
end

struct CooAllocation{T,A,B,C}
  counter::SparseMatrixCounter{T,A}
  I::B
  J::B
  V::C
end

LoopStyle(::Type{<:CooAllocation}) = Loop()

@inline function add_entry!(::typeof(+),a::CooAllocation{T},::Nothing,i,j) where T
  if is_entry_stored(T,i,j)
    a.counter.nnz = a.counter.nnz + 1
    k = a.counter.nnz
    a.I[k] = i
    a.J[k] = j
  end
  nothing
end

@inline function add_entry!(::typeof(+),a::CooAllocation{T},v,i,j) where T
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
  nz[k] = combine(v,Aij)
  A
end

function copy_entries!(a::T,b::T) where T<:AbstractSparseMatrix
  na = nonzeros(a)
  nb = nonzeros(b)
  if na !== nb
    copyto!(na,nb)
  end
end

function fill_entries!(A::AbstractSparseMatrix,v)
  nonzeros(A) .= v
  A
end

function allocate_coo_vectors(
   ::Type{<:AbstractSparseMatrix{Tv,Ti}},n::Integer) where {Tv,Ti}
  (zeros(Ti,n), zeros(Ti,n), zeros(Tv,n))
end

