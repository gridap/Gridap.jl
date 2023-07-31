# Block Assemblers

This file contains implementation details for block assemblers. We will first have a look at how standard sparse matrix assembly works, then compare it to block assembly.

## SparseMatrixAssemblers

Let's understand the main workflow for `SparseMatrixAssemblers` by looking at the `assemble_matrix` method:

```julia
  function assemble_matrix(a::SparseMatrixAssembler,matdata)
1   m1 = nz_counter(get_matrix_builder(a),(get_rows(a),get_cols(a)))
2   symbolic_loop_matrix!(m1,a,matdata)
3   m2 = nz_allocation(m1)
4   numeric_loop_matrix!(m2,a,matdata)
5   m3 = create_from_nz(m2)
    return m3
  end
```

By line number:

1) We retrieve the `SparseMatrixBuilder/ArrayBuilder` from the assembler. This structure has information on which type of array will be allocated at the end. For instance CSR vs CSC sparse matrix storage, which type of element type the array should hold (Float64, Float32,...), etc... With this information, we create a non-zero counter object `m1`, which will manage the counting of sparse entries and eventually the allocation of the array.
2) We do a symbolic loop over the matrix contributions in `matdata`. During this loop, we save the row/col indexes of the non-zeros that will be present in the final matrix.
3) We allocate the necessary space for the nonzero entries we have counted. In the case of CSR and CSC storage types, we do not allocate the final matrix yet but rather do everything in COO format (which is more efficient for random access assembly). `m2` is now an array allocator object.
4) We do a second loop over the matrix contributions in `matdata`, During this loop, we actually assemble the entries that will end up in the matrix.
5) We use the matrix allocator to allocate and fill the final sparse matrix.

So the objects involved and the overall workflow is given by

```md
SparseMatrixBuilder -> ArrayCounter -> ArrayAllocator -> Array
```

The second part of the puzzle is given by the loops over the data, i.e `symbolic_loop_X!` and `numeric_loop_X!`. Both loops are quite similar, so we will focus on the numeric loop, which is implemented in the following function:

```julia
  function numeric_loop_matrix!(A,a::SparseMatrixAssembler,matdata)
    get_mat(a::Tuple) = a[1]
    get_mat(a) = a
    if LoopStyle(A) == DoNotLoop()
      return A
    end
1   strategy = get_assembly_strategy(a)
    for (cellmat,_cellidsrows,_cellidscols) in zip(matdata...)
2     cellidsrows = map_cell_rows(strategy,_cellidsrows)
2     cellidscols = map_cell_cols(strategy,_cellidscols)
      @assert length(cellidscols) == length(cellidsrows)
      if length(cellidscols) > 0
        rows_cache = array_cache(cellidsrows)
        cols_cache = array_cache(cellidscols)
        mat1 = get_mat(first(cellmat))
        rows1 = getindex!(rows_cache,cellidsrows,1)
        cols1 = getindex!(cols_cache,cellidscols,1)
3       add! = AddEntriesMap(+)
3       add_cache = return_cache(add!,A,mat1,rows1,cols1)
3       caches = add_cache, vals_cache, rows_cache, cols_cache
3       _numeric_loop_matrix!(A,caches,cellmat,cellidsrows,cellidscols)
      end
    end
    A
  end

  @noinline function _symbolic_loop_matrix!(A,caches,cell_rows,cell_cols,mat1)
    touch_cache, rows_cache, cols_cache = caches
    touch! = TouchEntriesMap()
4   for cell in 1:length(cell_cols)
4     rows = getindex!(rows_cache,cell_rows,cell)
4     cols = getindex!(cols_cache,cell_cols,cell)
4     vals = getindex!(vals_cache,cell_vals,cell)
4     evaluate!(add_cache,add!,mat,vals,rows,cols)
4   end
  end
```

By line number:

1) We retrieve the `AssemblyStrategy` object from the assembler. This object contains all the information necessary to map the DoF ids from our mesh to the final columns/rows of the matrix. In serial this map is almost always the identity, but in parallel it is crucial to handle ghosts and local/global indexes.
2) We use the col/row maps in strategy to map the cell DoF ids in each cell contribution to the corresponding rows/cols in the final matrix.
3) We prepare a `TouchEntriesMap` (symbolic loop) or `AddEntriesMap` (numeric loop). These maps will be executed on each cell contribution, and are the ones responsible to allocate/assemble the contributions on the array counter.
4) For each cell, we retrieve the rows, cols and values and execute the map on these. This will allocate/assemble the contributions of this cell into the counter. In the case of `MultiFieldFESpaces`,  the cell indices `rows`/`cols` will be `VectorBlocks` with as many blocks as fields the `MultiFieldFESpace` has and `vals` will be a `MatrixBlock` with an array of blocks of size (# blocks in rows, # blocks in cols). The `TouchEntriesMap` and`AddEntriesMap` maps are specialised on these types, and assemble all blocks into the same `ArrayCounter`.

## BlockSparseMatrixAssemblers

To activate the block assemblers, we have created a new `MultiFieldStyle` called `BlockMultiFieldStyle`.  The purpose of this style is two-fold:

1) It activates the block assembly automatically when calling `SparseMatrixAssembler`, so that the everything fits with the high-level API.
2) It manages the numbering of the cell DoFs when performing the integrals, ensuring the DoFs ids are local to each block.

To create a block-assembled `MultiFieldFESpace`, you can use the following constructor:

```julia
mfs = BlockMultiFieldStyle()
Yb  = MultiFieldFESpace(tests;style=mfs)
Xb  = MultiFieldFESpace(trials;style=mfs)
```

By default, the final matrix and vector will have a block for each input `FESpace`. However, you can introduce some parameters when building your `BlockMultiFieldStyle` so that multiple fields are assembled in the same block (see section B).

### A) One to one Block <-> FESpace correspondence

The design of `BlockSparseMatrixAssemblers` is quite simple: Instead of having a single `SparseMatrixBuilder`, `ArrayCounter` and`Array` in which we assemble the entries coming from all fields in the `MultiFieldFESpace`, the assembler will create one of these objects for each final block and put them in a `ArrayBlock` object. We will then dispatch on the `ArrayBlock` type so that the contributions from each field is assembled in the block we want.

For instance, for the `nz_counter` function we dispatch as follows:

```julia
function Algebra.nz_counter(builder::MatrixBlock,axs)
  s = size(builder)
  rows = axs[1]
  cols = axs[2]
  counters = [nz_counter(builder.array[i,j],(rows[i],cols[j])) for i in 1:s[1], j in 1:s[2]]
  return ArrayBlock(counters,fill(true,size(counters)))
end
```

In this function, the variable `builder` is a `MatrixBlock{<:ArrayBuilder}`, which holds an array with the array builders for each final block.
We then simply select the rows/cols for each of the blocks and apply the `nz_counter` function to the corresponding `ArrayBuilder` ,  then return a `MatrixBlock{<:ArrayCounter}` which holds the array counters for each block.
Similar dispatches are provided for `nz_allocation`, `create_from_nz`, `map_cell_rows` and `map_cell_cols`.

We also specialise the evaluation of the `TouchEntriesMap` and`AddEntriesMap` maps when the counters are `BlockArrays`. For instance, let's have a look at the following function:

```julia
# A) Default implementation for MultiFieldFESpaces
function Fields.evaluate!(k::AddEntriesMap,A,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
  ni,nj = size(v.touched)
  for j in 1:nj
    for i in 1:ni
      if v.touched[i,j]
        evaluate!(cache[i,j],k,A,v.array[i,j],I.array[i],J.array[j])
      end
    end
  end
end

# B) Dispatching for block assemblers
function Fields.evaluate!(k::AddEntriesMap,A::MatrixBlock,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
  ni,nj = size(v.touched)
  for j in 1:nj
    for i in 1:ni
      if v.touched[i,j]
        evaluate!(cache[i,j],k,A.array[i,j],v.array[i,j],I.array[i],J.array[j])
      end
    end
  end
end
```

In the monolithic assembly of `MultiFieldFESpaces`, the variable `A` is an `ArrayCounter`. As you can see, all contribution blocks (coming from different fields) are assembled into the same `ArrayCounter`. The block-assembly counterpart will have the input `A` be a `MatrixBlock{<:ArrayCounter}`, and assembles each contribution block to it's corresponding `ArrayCounter` (notice the `A.array[i,j]`).

### B) Assembling multiple FE Fields into the same Block

The `BlockMultiFieldStyle` constructor can take up to three parameters:

1) `NB` :: Integer, representing the number of final blocks. Then the matrix and vector will have `NBxNB` and `NB` blocks respectively.
2) `SB` :: Tuple of integers, of length `NB`.  In each position, `SB[ib]` is the number of fields that will be assembled in that block.
3) `P`  :: Tuple of integers, of length the number of fields. This represents a field permutation, such that the fields will be reordered as`[P[1],P[2],....,P[n]]`.

Using this three parameters, one can assemble an arbitrary number of fields into any number of blocks.

**Example**: Consider we are solving an MHD problem with variables `(u,p,j,q)` , i.e (fluid velocity, fluid pressure, magnetic current, electric potential). Although the variables are in this specific order in the `MultiFieldFESpace`, we want to build a block-preconditioner that solves `(u,j)` together in a single block then `p` and `q` separately in two other blocks. Then we would need to assemble our system using `NB=3`, `SB=(2,1,1)` and `P=(1,3,2,4)`.
With this configuration, we will create 3 blocks. The first block will have size 2 and hold variables `[P[1],P[2]] = [1,3] = [u,j]`. The second block will have size 1 and hold variables `[P[3]] = [2] = [p]`. Finally, the third block will hold variables `[P[4]] = [4] = [q]`.

In terms of implementation, everything is the same. We use `ArrayBlockViews` (which is a view counterpart of `ArrayBlock`) so that an array of `NBxNB` array builders / array counters can be indexed using the field indexes. This allows us to use the same dispatches as we had in part A.
