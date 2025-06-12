
# Development notes for polytopal methods

## Hybrid assembly

The idea behind the following implementation is the following:

After a lot of thought, I think we probably cannot create cell-wise facet basis functions. The problem is efficiency: If we do so, we would be forced to pad the elemental cell matrices with zeros for those dofs evaluated at points that do not belong to the facet. We can still explore this idea for specific discretizations, but I think for methods where the dofs actually belong to the faces (i.e HDG, HHO, ...) this is not the way to go.

I think the key to efficiently implement hybrid assembly is leveraging the concept of patch-assembly, similar to what we do for GMG:

- First, we want to allow for overlapped triangulations. In this way, we can imagine a hybrid skeleton as a boundary triangulation of facets wrapped around the cells. I.e instead of what the current SkeletonTriangulation does, we have a single triangulation that contains each facet twice (once for each neighboring cell). The only thing missing is to tie back each facet to the cell it belongs to.
- This is where patches come into play. We have already been exploring patches for GMG + domain decomposition methods. In the simplest form, we consider one patch per cell, and each patch contains the cell itself and all touching facets.
- In the future, we could imagine more complex patches that contain more than one cell. This ties into the idea Santi had to explore small-domain BDDC solvers (where the coarse space is not tied to processors, but rather to some arbitrary domain decomposition).
- Once we have our patches, we have to create the machinery to assemble each patch individually. I.e we take the mesh-wide arrays, select the cells/faces that belong to each patch and assemble patch-local matrices. This can be done lazily (i.e we reuse the same matrix-vector for each patch).
- On top the of the resulting array of patch-local systems, one can easily create patch-local maps to implement features such as static condensations for HDG, reconstruction operators for HHO, etc...

The following are the files containing the implementation of the above ideas:

- Patches and patch triangulations: `src/Geometry/PatchTriangulations.jl`
- Patch Assembly for single-field: `src/FESpaces/PatchAssemblers.jl`
- Patch Assembly for multi-field: `src/MultiField/PatchAssemblers.jl`

## Logging some difficulties

For static condensation, the current approach is to gather the `cell_dof_ids` patch per patch, to create the `patch_dof_ids`. This gets rid of the face-wise splitting and merges all dofs belonging to a single patch into a single array.

Then, given a tuple `(dofs,patch)`, we can find the local ids of those dofs within the `patch_dofs[patch]` array. To do this efficiently, we have to sort each of the sub-arrays in `patch_dofs` in order to be able to use `searchsortedfirst`. Also, we remove completely boundary conditions, since the matrix-vector pairs resulting from static condensation already have had their boundary conditions applied.

The issue that this causes is that for HHO we will have two instances of the patch dofs: One with boundary conditions and one without, and both of them will have all dofs (i.e also with negative ids).
The one without boundary conditions is needed to create local projection operators. The one with is to globally assemble the result of the local projection operators.

So for a cell

```
        2
   ----------
   |        |
  3|        |-1
   |        |
   ----------
        1
```

we would have 3 versions of the patch dofs: 

```
  [1,2,3]    # bcs removed
  [1,2,3,4]  # space without bcs
  [1,2,3,-1] # space with bcs
```
