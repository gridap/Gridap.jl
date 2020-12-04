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
using Gridap.Algebra

using DocStringExtensions
using Test
using FillArrays
using LinearAlgebra
using BlockArrays
using Base: @propagate_inbounds
using ForwardDiff

import Base: size
import Base: getindex, setindex!
import Base: similar
import Base: IndexStyle

import Gridap.Algebra: scale_entries!
import Gridap.Algebra: fill_entries!

# CachedArray

export CachedArray
export CachedMatrix
export CachedVector
export setsize!
export setaxes!

# Map

export Map
export return_cache
export evaluate!
export evaluate
export test_map
export return_type
export return_value
export testargs
export inverse_map

export Broadcasting

export Operation


# LazyArray

export LazyArray
export array_cache
export getindex!
export lazy_map

# BlockArrays

export BlockArrayCoo
export BlockVectorCoo
export BlockMatrixCoo
export is_zero_block
export is_nonzero_block
export enumerateblocks
export eachblockid
export VectorOfBlockArrayCoo
export VectorOfBlockVectorCoo
export VectorOfBlockMatrixCoo
export zeros_like
export TwoLevelBlockedUnitRange
export append_ranges
export MultiLevelBlockedUnitRange
export blocks_equal
export num_blocks_equal
export local_range
export similar_range

export array_cache
export getindex!
export testitem
export testvalue
# export uses_hash
export test_array
export get_array
# # export add_to_array!


export CompressedArray

export Reindex
#export reindex

export PosNegReindex
export PosNegPartition
#export posneg_reindex

export FilterMap

export MulAddMap

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

export IdentityVector

export SubVector
export pair_arrays
export unpair_arrays

export AppendedArray
export lazy_append
export lazy_split

export autodiff_array_gradient
export autodiff_array_jacobian
export autodiff_array_hessian

export VectorWithEntryRemoved
export VectorWithEntryInserted

import AbstractTrees
export TreeNode
export print_op_tree
export similar_tree_node

import Gridap.Io: to_dict
import Gridap.Io: from_dict

include("Interface.jl")

include("Maps.jl")

include("AlgebraMaps.jl")

include("BlockArraysCoo.jl")

include("CachedArrays.jl")

include("LazyArrays.jl")

include("CompressedArrays.jl")

include("Tables.jl")

include("IdentityVectors.jl")

include("PosNegReindex.jl")

include("Reindex.jl")

include("FilteredArrays.jl")

include("SubVectors.jl")

include("ArrayPairs.jl")

include("AppendedArrays.jl")

include("Autodiff.jl")

include("VectorsWithEntryRemoved.jl")

include("VectorsWithEntryInserted.jl")

include("PrintOpTrees.jl")


end # module
