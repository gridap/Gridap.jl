
# FESpaceWithLinearConstraints

This is a space built from an initial FESpace plus a set of linear constraints.
The initial FESpace may also be a `FESpaceWithLinearConstraints`.

## Assumptions

- A constrained dof depends only on master dofs, i.e., only one level of constraints. This can be easily relaxed in the future.
- Free dofs can be constrained both by free and Dirichlet dofs. Dirichlet dofs can be constrained by Dirichlet dofs only.
- All constrained dofs belong to the original space. Master dofs will generally also belong to the original space, but they may be external. This allows us to have, e.g.
  - External master Dirichlet dofs, to implement affine constraints.
  - In distributed context, masters may belong to other processors and therefore NOT be local (and thus not belong to the original local space).

## Notation

We will setup some notation to differentiate between the different dof numberings that appear throughout the implementation:

We have two ways of numbering dofs:

- We will refer to `dof` in non-capital letters to refer to signed dof indices.
  The sign will indicate if the dof is free (positive) or Dirichlet (negative).
- We will refer to `DOF` in capital letters to refer to unsigned dof indices. In this case, Dirichlet dofs are also represented with a positive id "pas the end" of free dof ids. I.e., we have `DOF = dof` if `dof > 0` and `DOF = -dof + n_fdofs` if `dof < 0`, where `n_fdofs` is the number of free dofs in the original space.

We will differentiate between three different sets of dofs:

- First, DoFs in the original space will be denoted by `dof`/`DOF`. They can either be free (`fdof`) or Dirichlet (`ddof`).
- Master DoFs are the non-constrained degrees of freedom. They will be denoted by `mdof`/`mDOF`. They can either be free (`fmdof`) or Dirichlet (`dmdof`).
- Slave DoFs are the constrained degrees of freedom. They will be denoted by `sdof`/`sDOF`.

Moreover, we will sometimes need to refer to local dofs in a cell, for which we will use `ldof` for the original space and `lmdof` for the constrained space.
