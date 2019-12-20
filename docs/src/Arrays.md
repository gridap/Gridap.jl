
```@meta
CurrentModule = Gridap.Arrays
```
# Gridap.Arrays

```@docs
Arrays
```

## Extended AbstractArray interface

When implementing new array types, it can be needed some scratch data (e.g., allocating the output), when recovering an item from an array (typically if the array elements are non-isbits objects). To circumvent this, the user could provide the scratch data needed when getting an item. However, the Julia array interface does not support this approach. When calling `a[i]`, in order to get the element with index `i` in array `a`, there is no extra argument for the scratch data. In order to solve this problem, we add new methods to the `AbstractArray` interface of Julia. We provide default implementations to the new methods, so that any `AbstractArray` can be used with the extended interface. New array implementations can overload these default implementations to improve performance. The most important among the new methods is [`getindex!`](@ref), which allows to recover an item in the array by passing some scratch data.

The new methods are:
- [`getindex!(cache,a::AbstractArray,i...)`](@ref)
- [`array_cache(a::AbstractArray)`](@ref)
- [`uses_hash(::Type{<:AbstractArray})`](@ref)
- [`testitem(a::AbstractArray)`](@ref)

These methods can be stressed with the following function
- [`test_array`](@ref)

```@docs
getindex!(cache,a::AbstractArray,i...)
array_cache(a::AbstractArray)
uses_hash(::Type{<:AbstractArray})
testitem(a::AbstractArray)
test_array
```

## Working with several arrays at once

```@docs
getitems!
array_caches
testitems
```

## Creting lazy operation trees

```@docs
apply(f,a::AbstractArray...)
apply(::Type{T},f,a::AbstractArray...) where T
apply(f::AbstractArray,a::AbstractArray...)
apply(::Type{T},f::AbstractArray,a::AbstractArray...) where T
apply_all
```

### Operation kernels

```@docs
Kernel
apply_kernel!(cache,f,x...)
kernel_cache(f,x...)
kernel_return_type(f,x...)
test_kernel
```

### Other functions using kernels

```@docs
apply_kernel
apply_kernels!
kernel_caches
kernel_return_types
```

### Build-in kernels

```@docs
bcast
elem
contract
```
## Helper functions

```@docs
collect1d
```

## Concrete array implementations

### CachedArray

```@docs
CachedArray
CachedArray(a::AbstractArray)
CachedArray(T,N)
setsize!
CachedMatrix
CachedVector
```
### CompressedArray

```@docs
CompressedArray
CompressedArray(::AbstractArray,::AbstractArray)
```

### Table

```@docs
Table
Table(data::AbstractVector,ptrs::AbstractVector)
Table(a::AbstractVector{<:AbstractVector})
identity_table(::Type{T},::Type{P},l::Integer) where {T,P}
empty_table(::Type{T},::Type{P}, l::Integer) where {T,P}
get_ptrs_eltype(::Table{T,P}) where {T,P}
get_data_eltype(::Table{T,P}) where {T,P}
append_tables_globally(tables::Table{T,P}...) where {T,P}
append_tables_locally(offsets::NTuple, tables::NTuple)
append_tables_locally(tables::Table...)
rewind_ptrs!(ptrs::AbstractVector{<:Integer})
generate_data_and_ptrs(vv::AbstractVector{<:AbstractVector{T}}) where T
length_to_ptrs!(ptrs::AbstractArray{<:Integer})
append_ptrs(pa::AbstractVector{T},pb::AbstractVector{T}) where T
append_ptrs!(pa::AbstractVector{T},pb::AbstractVector{T}) where T
find_inverse_index_map(a_to_b, nb)
find_inverse_index_map!(a_to_b, b_to_a)
flatten_partition(a_to_bs::Table,nb::Integer)
find_local_index(a_to_b, b_to_la_to_a)
get_local_item(a_to_lb_to_b, lb::Integer)
UNSET
```
### LocalToGlobalArray

```@docs
LocalToGlobalArray
LocalToGlobalArray(::AbstractArray{<:AbstractArray},::AbstractArray)
```

