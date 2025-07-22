
# Gridap ecosystem

The Gridap ecosystem has grown over the years, and now includes many packages that vastly extend the capabilities of the base library.
We here list the ones that are currently maintained by the Gridap team, as well as external packages maintained by the community.

!!! note

    We would also love to showcase any projects (registered or not) that use Gridap. **If you are developper of a project and wish to be listed here, please contact us and we'll be happy to add your project to the list!**

## Meshing

- [`GridapGmsh.jl`](https://github.com/gridap/GridapGmsh.jl) is our interface to the [Gmsh](https://gmsh.info/) mesh generator. It allows you to create meshes from julia scripts, or import existing meshes in the Gmsh format.

## Visualization

- [`GridapMakie.jl`](https://github.com/gridap/GridapMakie.jl) is our interface to the [Makie](https://makie.juliaplots.org/stable/) visualization library.

## Distributed computing

Gridap supports distributed-memory computing through MPI:

- [PartitionedArrays.jl](https://github.com/PartitionedArrays/PartitionedArrays.jl) provides distributed algebra interfaces.
- [GridapDistributed.jl](https://github.com/gridap/GridapDistributed.jl) extends the Gridap API for distributed computing. This package preserves (almost perfectly) the high-level interface of Gridap, allowing user to almost seamlessly parallelize their code accross multiple machines.

## Solvers

Beyond the basic solvers provided by the base library, Gridap has several packages providing advanced solvers for distributed computing:

- [GridapSolvers.jl](https://github.com/gridap/GridapSolvers.jl) provides physics-informed solvers and preconditioners for HPC, built fully in Julia. Among other things, it includes:
  - A collection of Krylov solvers
  - An interface to design block preconditioners
  - Geometric Multigrid solvers
- [`GridapPETSc.jl`](https://github.com/gridap/GridapPETSc.jl) is our interface to the [PETSc](https://petsc.org/release/) library, which provides a wide range of solvers and preconditioners. These can also be used in conjunction with the block-preconditioners of `GridapSolvers.jl`.
- [`GridapPardiso.jl`](https://github.com/gridap/GridapPardiso.jl) is our interface to the Pardiso multithreaded solver. A newer interface based on the official package [`Pardiso.jl`](https://github.com/JuliaSparse/Pardiso.jl) is also available as an extension of `GridapSolvers.jl`.

## Embedded methods

Embedded methods are a popular way of dealing with complex geometries. We provide a high order interface for these methods in the [`GridapEmbedded.jl`](https://github.com/gridap/GridapEmbedded.jl) package. This package allows you to define geometries using functions and level-sets. Support for STL meshes is provides in the satellite package [`STLCutters.jl`](https://github.com/gridap/STLCutters.jl).

## Adaptive mesh refinement

Gridap provides a high-level expandable API for refining and coarsening meshes, as well as support for AMR.
We also provide support for the well-known library [`p4est`](https://www.p4est.org/) through the package [`GridapP4est.jl`](https://github.com/gridap/GridapP4est.jl), which provides non-conforming highly efficient AMR on parallel distributed memory machines.

## Topology optimization

- [GridapTopOpt.jl](https://github.com/zjwegert/GridapTopOpt.jl) is a package for topology optimization using level-set methods.

## Reduced order methods

- [GridapROM.jl](https://github.com/gridap/GridapROMs.jl) provides tools for the solution of parameterized partial differential equations (PDEs) with reduced order models (ROMs).

## Applications

Some repositories providing applications using `Gridap`:

- [GridapGeosciences.jl](https://github.com/gridapapps/GridapGeosciences.jl), a package to solve geophysical flow problems.
- [GridapMHD.jl](https://github.com/gridapapps/GridapMHD.jl), a package for HPC-focussed magneto-hydrodynamics simulations.
