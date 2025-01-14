
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
