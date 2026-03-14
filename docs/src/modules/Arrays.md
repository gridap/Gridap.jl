# Gridap.Arrays

The `Gridap.Arrays` module provides a high-performance array and mapping infrastructure designed to minimize allocations through lazy evaluation and mutable caches. This module is the engine behind many of Gridap's efficient operations.

The main abstractions are the `Map` interface (functional engine) and the extended `AbstractArray` interface (efficient iteration). Key data structures include:

- `LazyArray`: Represents the result of an operation on other arrays, computed only on demand.
- `CompressedArray`: Memory-efficient storage for arrays with many repeated values.
- `Fill`: (from `FillArrays.jl`) Efficient representation of constant data across an entire array.
- `Table`: A memory-efficient representation of jagged arrays (e.g., cell-to-node connectivity).

```@docs
Arrays
```

#### Contents

```@contents
Pages = ["Arrays.md"]
Depth = 2:3
```

## Map Interface

The `Map` interface is the functional engine of Gridap. It defines objects that provide a cache and an in-place evaluation for performance, specially useful in lazy operations.

```@docs
Map
return_cache
evaluate!
evaluate
test_map
return_type
```

### Operations and broadcasting

The following types represent common higher-level operations on maps.

```@docs
Broadcasting
Operation
```

### Reindex maps

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/Reindex.jl","/PosNegReindex.jl"]
```

### Algebra maps

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/AlgebraMaps.jl"]
```

### Other maps

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/KeyToValMaps.jl"]
```

## Array Interface

Gridap extends the `AbstractArray` interface to support efficient iteration through mutable caches.

```@docs
array_cache
getindex!
testitem
test_array
invalidate_cache!
```

### LazyArrays

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/LazyArrays.jl"]
```

### CompressedArrays

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/CompressedArrays.jl"]
```

### Tables

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/Tables.jl"]
```

### CachedArrays

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/CachedArrays.jl"]
```

### Other arrays

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/IdentityVectors.jl","/FilteredArrays.jl","/AppendedArrays.jl","/VectorsWithEntryRemoved.jl","/VectorsWithEntryInserted.jl","/ArrayPairs.jl"]
```

## ArrayBlocks

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/ArrayBlocks.jl"]
```

## Automatic Differentiation

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/Autodiff.jl"]
```

## Printing operation trees

```@autodocs
Modules = [Arrays,]
Order   = [:type, :constant, :macro, :function]
Pages   = ["/PrintOpTrees.jl"]
```
