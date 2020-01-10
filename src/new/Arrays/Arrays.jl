"""
This module provides:
- An extension of the `AbstractArray` interface in order to properly deal with mutable caches.
- A mechanism to generate lazy arrays resulting from operations between arrays.
- A collection of concrete implementations of `AbstractArray`.

The exported names in this module are:

$(EXPORTS)
"""
module Arrays

using Gridap.Helpers
using Gridap.Inference

using DocStringExtensions
using Test
using FillArrays
using Base: @propagate_inbounds

export array_cache
export getindex!
export getitems!
export testitem
export uses_hash
export test_array
export testitems
export array_caches
export get_array

export CachedArray
export CachedMatrix
export CachedVector
export setsize!

export CompressedArray
export LocalToGlobalArray
export LocalToGlobalPosNegArray

export kernel_cache
export kernel_caches
export apply_kernels!
export apply_kernel!
export apply_kernel
export test_kernel
export bcast
export elem
export contract
export kernel_return_type
export kernel_return_types
export Kernel

export apply
export apply_all

export Table
export identity_table
export empty_table
export rewind_ptrs!
export length_to_ptrs!
export append_ptrs
export append_ptrs!
export get_ptrs_eltype
export get_data_eltype
export generate_data_and_ptrs
export find_inverse_index_map
export find_inverse_index_map!
export append_tables_globally
export append_tables_locally
export flatten_partition
export collect1d
export UNSET
export get_local_item
export find_local_index

export reindex

import Base: size
import Base: getindex, setindex!
import Base: similar
import Base: IndexStyle

import Gridap.Io: to_dict
import Gridap.Io: from_dict

include("Interface.jl")

include("CachedArrays.jl")

include("Kernels.jl")

include("Apply.jl")

include("CompressedArrays.jl")

include("Tables.jl")

include("LocalToGlobalArrays.jl")

include("LocalToGlobalPosNegArrays.jl")

include("Reindex.jl")

end # module
