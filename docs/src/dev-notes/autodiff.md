
# Autodiff

## How it works

Gridap makes use of [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/) to take derivatives of FE-based functionals. The main idea is to propagate dual numbers through the FE assembly process, then extract the partials from the dualized result. For a more detailed explanation, please check out `ForwardDiff.jl`'s documentation.

Due to the cell-wise (local) nature of FE integral operators, we can compute the contributions of each cell to the gradient/jacobian independently, then assemble the local contributions into the global gradient/jacobian. This means we can re-use all of the existing FE evaluation and assembly machinery. This also means we can use the number of DoFs per cell as a very natural choice for our chunk size.

Let's consider a simple scalar-valued functional ``j(u)`` that takes as argument an `FEFunction` ``u`` defined on an `FESpace` ``U``. We want to compute the gradient of ``j`` w.r.t ``u`` at a certain point ``u_0``, i.e ``\nabla j(u_0)``. This is a vector of size `n = num_free_dofs(U)`, where each entry is the partial derivative of ``j`` w.r.t the corresponding DoF of ``u``, i.e ``\nabla j(u_0)|_i = (\partial j / \partial u_i)(u_0)``.

Evaluating this gradient involves the following steps:

- First, we seed the dual numbers as cell-wise dof values for ``u_0``, i.e in each cell we set ``u^* = \sum dual(u_i) \phi_i`` where ``phi_i`` are the basis functions of `U` within teh cell. Note that the seeds are generated independently in each cell, i.e ``dual(u_i)`` will have as many partials as the number of DoFs in the cell.
- We then evaluate the functional ``j(u^*)`` in each cell, but do not assemble the result. This gives us a cell-wise array of dual numbers.
- We extract the partials in each cell, which produces a cell-wise array of vectors. The local vectors hold the partials of ``j`` w.r.t the DoFs of ``u`` in the cell. We then assemble these local vectors into the global gradient, using the FE assembly machinery.

Note this process is fully local, i.e we can do all of the above in a single cell, then move to the next cell. We can also use the lazy-evaluation machinery of `Gridap` to reuse all the evaluation caches.

## Implementation

The above process is implemented in the following function (from `src/Arrays/Autodiff.jl`):

```julia
function autodiff_array_gradient(a,i_to_x)
  dummy_tag = ()->()
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient,dummy_tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x) # Dualize cell dofs
  i_to_ydual = a(i_to_xdual) # Evaluate the functional
  i_to_result = lazy_map(AutoDiffMap(),i_to_cfg,i_to_ydual) # Extract partials
  i_to_result
end
```

where `i_to_x` has the cell dof values of the evaluation point ``u_0`` and `a` is an auxiliary function that takes an array of dualized cell dofs and evaluates the functional provided by the user, i.e (for the above example):

```julia
a(dual_cell_dofs) = j(FEFunction(U,dual_cell_dofs))
```

The intermediate array `i_to_cfg` holds the `ForwardDiff` seeding configurations for each cell. The resulting array `i_to_result` is then assembled into the global gradient.

Sometimes, the functional ``j`` is integrated on a `Triangulation` different from the one where ``U`` is defined. This case requires an extra step, as follows:

```julia
function autodiff_array_gradient(a,i_to_x,j_to_i)
  dummy_tag = ()->()
  i_to_cfg = lazy_map(ConfigMap(ForwardDiff.gradient,dummy_tag),i_to_x)
  i_to_xdual = lazy_map(DualizeMap(),i_to_cfg,i_to_x)
  j_to_ydual = a(i_to_xdual)
  j_to_cfg = autodiff_array_reindex(i_to_cfg,j_to_i)
  j_to_result = lazy_map(AutoDiffMap(),j_to_cfg,j_to_ydual)
  j_to_result
end
```

where `j_to_i` is a mapping from the cells of the integrand `Triangulation` to the cells of the `FESpace` `U`. We then move the configurations into the final `Triangulation`, to be able to correctly extract the partials. The remaining steps are the same as before.

Analog functions are implemented for the Jacobian, where local contributions will now be matrices that hold ``(\partial r_j / \partial u_i)(u_0)`` for each local residual ``r_j`` (corresponding to the local DoFs of the test space) and each local DoF of the trial space. These contributions are then assembled into the global Jacobian.
