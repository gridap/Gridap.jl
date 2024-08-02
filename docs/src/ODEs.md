```@meta
CurrentModule = Gridap.ODEs
```

# Gridap.ODEs

We consider an initial value problem written in the form
```math
\left\{\begin{array}{rcll}\boldsymbol{r}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n} \boldsymbol{u}) &=& \boldsymbol{0}_{d}, \\ \partial_{t}^{k} \boldsymbol{u}(t_{0}) &=& \boldsymbol{u}_{0}^{k} & 0 \leq k \leq n-1,\end{array}\right.
```

where

* ``\boldsymbol{u}: \mathbb{R} \to \mathbb{R}^{d}`` is the unknown of the problem,
* ``n \in \mathbb{N}`` is the order of the ODE,
* ``t_{0} \in \mathbb{R}`` is the initial time and ``\{\boldsymbol{u}_{0}^{k}\}_{0 \leq k \leq n-1} \in (\mathbb{R}^{d})^{n-1}`` are the initial conditions, and
* ``\boldsymbol{r}: \mathbb{R} \times (\mathbb{R}^{d})^{n} \to \mathbb{R}^{d}`` is the residual.

> We illustrate these notations on the semi-discretisation of the heat equation. It is a first-order ODE so we have ``n = 1``, and ``d`` is the number of degrees of freedom. The residual and initial condition have the form
> ```math
> \boldsymbol{r}(t, \boldsymbol{u}, \dot{\boldsymbol{u}}) \doteq \boldsymbol{M} \dot{\boldsymbol{u}} + \boldsymbol{K}(t) \boldsymbol{u} - \boldsymbol{f}(t), \qquad \boldsymbol{u}(t_{0}) = \boldsymbol{u}_{0}^{0},
> ```
> where ``\boldsymbol{M} \in \mathbb{R}^{d \times d}`` is the mass matrix, ``\boldsymbol{K}: \mathbb{R} \to \mathbb{R}^{d \times d}`` is the stiffness matrix, ``\boldsymbol{f}: \mathbb{R} \to \mathbb{R}^{d}`` is the forcing term, and ``\boldsymbol{u}_{0}^{0}`` is the initial condition.

Suppose that we are willing to approximate ``\boldsymbol{u}`` at a time ``t_{F} > t_{0}``. A numerical scheme splits the time interval ``[t_{0}, t_{F}]`` into smaller intervals ``[t_{n}, t_{n+1}]`` (that do not have to be of equal length) and propagates the information at time ``t_{n}`` to time ``t_{n+1}``. More formally, we consider a general framework consisting of a starting, an update, and a finishing map defined as follows
* The **starting map** ``\mathcal{I}: (\mathbb{R}^{d})^{n} \to (\mathbb{R}^{d})^{s}`` converts the initial conditions into ``s`` state vectors, where ``s \geq n``.
* The **marching map** ``\mathcal{U}: \mathbb{R} \times (\mathbb{R}^{d})^{s} \to (\mathbb{R}^{d})^{s}`` updates the state vectors from time ``t_{n}`` to time ``t_{n+1}``.
* The **finishing map** ``\mathcal{F}: \mathbb{R} \times (\mathbb{R}^{d})^{s} \to \mathbb{R}^{d}`` converts the state vectors into the evaluation of ``\boldsymbol{u}`` at the current time.

> In the simplest case, the time step ``h = h_{n} = t_{n+1} - t_{n}`` is prescribed and constant across all iterations. The state vectors are simply the initial conditions, i.e. ``s = n`` and ``\mathcal{I} = \mathrm{id}``, and assuming that the initial conditions are given by increasing order of time derivative, ``\mathcal{F}`` returns its first input.
>
> Some schemes need nontrivial starting and finishing maps. (See the generalised-``\alpha`` schemes below.) When higher-order derivatives can be retrieved from the state vectors, it is also possible to take another definition for ``\mathcal{F}`` so that it returns the evaluation of ``\boldsymbol{u}`` and higher-order derivatives at the current time.

These three maps need to be designed such that the following recurrence produces approximations of ``\boldsymbol{u}`` at the times of interest ``t_{n}``
```math
\left\{\begin{array}{lcl}
\left\{\boldsymbol{s}\right\}_{n+1} &=& \mathcal{U}(h_{n}, \left\{\boldsymbol{s}\right\}_{n}) \\
\left\{\boldsymbol{s}\right\}_{0} &=& \mathcal{I}(\boldsymbol{u}_{0}^{0}, \ldots, \boldsymbol{u}_{0}^{n-1})
\end{array}\right., \qquad \boldsymbol{u}_{n} = \mathcal{F}(\left\{\boldsymbol{s}\right\}_{n}).
```
More precisely, we would like ``\boldsymbol{u}_{n}`` to be close to ``\boldsymbol{u}(t_{n})``. Here the notation ``\{\boldsymbol{s}\}_{n}`` stands for the state vector, i.e. a vector of ``s`` vectors: ``\{\boldsymbol{s}\}_{n} = (\boldsymbol{s}_{n, i})_{1 \leq i \leq s}``. In particular, we notice that we need the exactness condition
```math
\mathcal{F} \circ \mathcal{I}(\boldsymbol{u}_{0}^{0}, \ldots, \boldsymbol{u}_{0}^{n-1}) = \boldsymbol{u}_{0}.
```
This is a condition on the design of the pair (``\mathcal{I}``, ``\mathcal{F}``).

# Classification of ODEs and numerical schemes
Essentially, a numerical scheme converts a (continuous) ODE into (discrete) nonlinear systems of equations. These systems of equations can be linear under special conditions on the nature of the ODE and the numerical scheme. Since numerical methods for linear and nonlinear systems of equations can be quite different in terms of cost and implementation, we are interested in solving linear systems whenever possible. This leads us to perform the following classifications.

## Classification of ODEs
We define a few nonlinearity types based on the expression of the residual.
* **Nonlinear**. Nothing special can be said about the residual.
* **Quasilinear**. The residual is linear with respect to the highest-order time derivative and the corresponding linear form may depend on time and lower-order time derivatives, i.e.
```math
\boldsymbol{r}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n} \boldsymbol{u}) = \boldsymbol{M}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n-1} \boldsymbol{u}) \partial_{t}^{n} \boldsymbol{u} + \boldsymbol{f}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n-1} \boldsymbol{u}).
```
We call the matrix ``\boldsymbol{M}: \mathbb{R} \to \mathbb{R}^{d \times d}`` the mass matrix. In particular, a quasilinear ODE is a nonlinear ODE.
* **Semilinear**. The residual is quasilinear and the mass matrix may only depend on time, i.e.
```math
\boldsymbol{r}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n} \boldsymbol{u}) = \boldsymbol{M}(t) \partial_{t}^{n} \boldsymbol{u} + \boldsymbol{f}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n-1} \boldsymbol{u}).
```
In particular, a semilinear ODE is a quasilinear ODE.
* **Linear**. The residual is linear with respect to all time derivatives, i.e.
```math
\boldsymbol{r}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n} \boldsymbol{u}) = \sum_{0 \leq k \leq n} \boldsymbol{A}_{k}(t) \partial_{t}^{k} \boldsymbol{u} - \boldsymbol{f}(t).
```
We refer to the matrix ``\boldsymbol{A}_{k}: \mathbb{R} \to \mathbb{R}^{d \times d}`` as the ``k``-th linear form of the residual. We may still define the mass matrix ``\boldsymbol{M} = \boldsymbol{A}_{n}``. Note that the term independent of $u$, i.e. the forcing term, is subtracted from the residual. This aligns with standard conventions, and in particular with those of `AffineFEOperator` (see example in _Finite element operators_ below, in the construction of a `TransientLinearFEOperator`). In particular, a linear ODE is a semilinear ODE.

> Note that for residuals of order zero (i.e. "standard" systems of equations), the definitions of quasilinear, semilinear, and linear coincide.

We consider an extra ODE type that is motivated by stiff problems. We say that an ODE has an implicit-explicit (IMEX) decomposition if it be can written as the sum of a residual of order ``n`` and another residual of order ``n-1``, i.e.
```math
\boldsymbol{r}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n} \boldsymbol{u}) = \boldsymbol{r}_{\text{implicit}}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n} \boldsymbol{u}) + \boldsymbol{r}_{\text{explicit}}(t, \partial_{t}^{0} \boldsymbol{u}, \ldots, \partial_{t}^{n-1} \boldsymbol{u}).
```
The decomposition takes the form above so that the mass matrix of the global residual is fully contained in the implicit part. The table below indicates the type of the corresponding global ODE.

| Explicit \ Implicit | Nonlinear   | Quasilinear | Semilinear | Linear     |
|---------------------|-------------|-------------|------------|------------|
| Nonlinear           | Nonlinear   | Quasilinear | Semilinear | Semilinear |
| Linear              | Nonlinear   | Quasilinear | Semilinear | Linear     |

In particular, for the global residual to be linear, both the implicit and explicit parts need to be linear too.

> In the special case where the implicit part is linear and the explicit part is quasilinear or semilinear, we could, in theory, identify two linear forms for the global residual. However, introducing this difference would call for an order-dependent classification of ODEs and this would create (infinitely) many new types. Since numerical schemes can rarely take advantage of this extra structure in practice, we still say that the global residual is semilinear in these cases.

## Classification of numerical schemes
We introduce a classification of numerical schemes based on where they evaluate the residual during the state update.

* If it is possible (up to a change of variables) to write the system of equations for the state update as evaluations of the residual at known values (that depend on the solution at the current time) for all but the highest-order derivative, we say that the scheme is explicit.
* Otherwise, we say that the scheme is implicit.

> For example, when solving a first-order ODE, the state update would involve solving one or more equations of the type 
> ```math
> \boldsymbol{r}(t_{k}, \boldsymbol{u}_{k}(\boldsymbol{x}), \boldsymbol{v}_{k}(\boldsymbol{x})) = \boldsymbol{0},
> ```
> where ``\boldsymbol{x}`` and the unknown of the state update. The scheme is explicit if it is possible to introduce a change of variables such that ``\boldsymbol{u}_{k}`` does not depend on ``\boldsymbol{x}``. Otherwise, it is implicit.

## Classification of systems of equations
It is advantageous to introduce this classification of ODE and numerical schemes because the system of equations arising from the discretisation of the ODE by a numerical scheme will be linear or nonlinear depending on whether the scheme is explicit, implicit, or implicit-explicit, and on the type of the ODE. More precisely, we have the following table.

|                   | Nonlinear   | Quasilinear | Semilinear | Linear |
|-------------------|-------------|-------------|------------|--------|
| Explicit          | Nonlinear   | Linear      | Linear     | Linear |
| Implicit          | Nonlinear   | Nonlinear   | Nonlinear  | Linear |

When the system is linear, another important practical consideration is whether the matrix of the system is constant across iterations or not. This is important because a linear solver typically performs a factorisation of the matrix, and this operation may only be performed once if the matrix is constant.
* If the linear system comes from an explicit scheme, the matrix of the system is constant if the mass matrix is. This means that the ODE has to be quasilinear.
* If the linear system comes from an implicit scheme, all the linear forms must be constant for the system to have a constant matrix.

## Reuse across iterations
For performance reasons, it is thus important that the ODE be described in the most specific way. In particular, we consider that the mass term of a quasilinear ODE is not constant, because if it is, the ODE is semilinear. We enable the user to specify the following constant annotations:
* For nonlinear and quasilinear ODE, no quantity can be described as constant.
* For a semilinear ODE, whether the mass term is constant.
* For a linear ODE, whether all the linear forms are constant.

If a linear form is constant, regardless of whether the numerical scheme relies on a linear or nonlinear system, it is always possible to compute the jacobian of the residual with respect to the corresponding time derivative only once and retrieve it in subsequent computations of the jacobian.

# High-level API in Gridap
The ODE module of `Gridap` relies on the following structure.

## Finite element spaces
The time-dependent counterpart of `TrialFESpace` is `TransientTrialFESpace`. It is built from a standard `TestFESpace` and is equipped with time-dependent Dirichlet boundary conditions.
> By definition, test spaces have zero Dirichlet boundary conditions so they need not be seen as time-dependent objects.

A `TransientTrialFESpace` can be evaluated at any time derivative order, and the corresponding Dirichlet values are the time derivatives of the Dirichlet boundary conditions.

For example, the following creates a transient `FESpace` and evaluates its first two time derivatives.
```
g(t) = x -> x[1] + x[2] * t
V = FESpace(model, reffe, dirichlet_tags="boundary")
U = TransientTrialFESpace (V, g)

t0 = 0.0
U0 = U(t0)

∂tU = ∂t(U)
∂tU0 = ∂tU(t0)

∂ttU = ∂tt(U) # or ∂ttU = ∂t(∂t(U))
∂ttU0 = ∂ttU(t0)
```

## Cell fields
The time-dependent equivalent of `CellField` is `TransientCellField`. It stores the cell field itself together with its derivatives up to the order of the ODE.

For example, the following creates a `TransientCellField` with two time derivatives.
```
u0 = zero(get_free_dof_values(U0))
∂tu0 = zero(get_free_dof_values(∂tU0))
∂ttu0 = zero(get_free_dof_values(∂ttU0))
u = TransientCellField(u0, (∂tu0, ∂ttu0))
```

## Finite element operators
The time-dependent analog of `FEOperator` is `TransientFEOperator`. It has the following constructors based on the nonlinearity type of the underlying ODE.

* `TransientFEOperator(res, jacs, trial, test)` and `TransientFEOperator(res, trial, test; order)` for the version with automatic jacobians. The residual is expected to have the signature `residual(t, u, v)`.
* `TransientQuasilinearFEOperator(mass, res, jacs, trial, test)` and `TransientQuasilinearFEOperator(mass, res, trial, test; order)` for the version with automatic jacobians. The mass and residual are expected to have the signatures `mass(t, u, dtNu, v)` and `residual(t, u, v)`, i.e. the mass is written as a linear form of the highest-order time derivative `dtNu`. In this setting, the mass matrix is supposed to depend on lower-order time derivatives, so `u` is provided for the nonlinearity of the mass matrix.
* `TransientSemilinearFEOperator(mass, res, jacs, trial, test; constant_mass)` and `TransientSemilinearFEOperator(mass, res, trial, test; order, constant_mass)` for the version with automatic jacobians. (The jacobian with respect to ``\partial_{t}^{n} \boldsymbol{u}`` is simply the mass term). The mass and residual are expected to have the signatures `mass(t, dtNu, v)` and `residual(t, u, v)`, where here again `dtNu` is the highest-order derivative. In particular, the mass is specified as a linear form of `dtNu`.
* `TransientLinearFEOperator(forms, res, jacs, trial, test; constant_forms)` and `TransientLinearFEOperator(forms, res, trial, test; constant_forms)` for the version with automatic jacobians. (In fact, the jacobians are simply the forms themselves). The forms and residual are expected to have the signatures `form_k(t, dtku, v)` and `residual(t, v)`, i.e. `form_k` is a linear form of the ``k``-th order derivative, and the residual does not depend on `u`.

It is important to note that all the terms are gathered in the residual, including the forcing term. In the common case where the ODE is linear, the residual is only the forcing term, and it is subtracted from the bilinear forms (see example below).

Here, in the signature of the residual, `t` is the time at which the residual is evaluated, `u` is a function in the trial space, and `v` is a test function. Time derivatives of `u` can be included in the residual via the `∂t` operator, applied as many times as needed, or using the shortcut `∂t(u, N)`.

Let us take the heat equation as an example. The original ODE is
```math
\partial_{t} u - \nabla \cdot (\kappa(t) \nabla u) = f(t),
```
where ``\kappa`` is the (time-dependent) thermal conductivity and ``f`` is the forcing term. We readily obtain the weak form
```math
\int_{\Omega} v \partial_{t} u(t) + \nabla v \cdot (\kappa(t) \nabla u(t)) \ \mathrm{d} \Omega = \int_{\Omega} v f(t) \ \mathrm{d} \Omega.
```
It could be described as follows.

* As a `TransientFEOperator`:
```
res(t, u, v) = ∫( v ⋅ ∂t(u) + ∇(v) ⋅ (κ(t) ⋅ ∇(u)) - v ⋅ f(t) ) dΩ
TransientFEOperator(res, U, V)
```
* As a `TransientQuasilinearFEOperator`:
```
mass(t, u, dtNu, v) = ∫( v ⋅ dtNu ) dΩ
res(t, u, v) = ∫( ∇(v) ⋅ (κ(t) ⋅ ∇(u)) - v ⋅ f(t) ) dΩ
TransientQuasilinearFEOperator(mass, res, U, V)
```
* As a `TransientSemilinearFEOperator`:
```
mass(t, dtu, v) = ∫( v ⋅ dtu ) dΩ
res(t, u, v) = ∫( ∇(v) ⋅ (κ(t) ⋅ ∇(u)) - v ⋅ f(t) ) dΩ
TransientSemilinearFEOperator(mass, res, U, V, constant_mass=true)
```
* As a `TransientLinearFEOperator`:
```
stiffness(t, u, v) = ∫( ∇(v) ⋅ (κ(t) ⋅ ∇(u)) ) dΩ
mass(t, dtu, v) = ∫( v ⋅ dtu ) dΩ
res(t, u, v) = ∫( v ⋅ f(t) ) dΩ
TransientLinearFEOperator((stiffness, mass), res, U, V, constant_forms=(false, true))
```
If ``\kappa`` is constant, the keyword `constant_forms` could be replaced by `(true, true)`.

## The `TimeSpaceFunction` constructor
Apply differential operators on a function that depends on time and space is somewhat cumbersome. Let `f` be a function of time and space, and `g(t) = x -> f(t, x)` (as in the prescription of the boundary conditions `g` above). Applying the operator ``\partial_{t} - \Delta``  to `g` and evaluating at ``(t, x)`` is written `∂t(g)(t)(x) - Δ(g(t))(x)`.

The constructor `TimeSpaceFunction` allows for simpler notations: let `h = TimeSpaceFunction(g)`. The object `h` is a functor that supports the notations 
* `op(h)`: a `TimeSpaceFunction` representing both `t -> x -> op(f)(t, x)` and `(t, x) -> op(f)(t, x)`,
* `op(h)(t)`: a function of space representing `x -> op(f)(t, x)`
* `op(h)(t, x)`: the quantity `op(f)(t, x)` (this notation is equivalent to `op(h)(t)(x)`),

for all spatial and temporal differential operator, i.e. `op` in `(time_derivative, gradient, symmetric_gradient, divergence, curl, laplacian)` and their symbolic aliases (`∂t`, `∂tt`, `∇`, ...). The operator above applied to `h` and evaluated at `(t, x)` can be conveniently written `∂t(h)(t, x) - Δ(h)(t, x)`.

## Solver and solution
The next step is to choose an ODE solver (see below for a full list) and specify the boundary conditions. The solution can then be iterated over until the final time is reached.

For example, to use the ``\theta``-method with a nonlinear solver, one could write
```
t0 = 0.0
tF = 1.0
dt = 0.1
uh0 = interpolate_everywhere(t0, U(t0))

res(t, u, v) = ∫( v ⋅ ∂t(u) + ∇(v) ⋅ (κ(t) ⋅ ∇(u)) - v ⋅ f(t) ) dΩ
jac(t, u, du, v) = ∫( ∇(v) ⋅ (κ(t) ⋅ ∇(du)) ) dΩ
jac_t(t, u, dtu, v) = ∫( v ⋅ dtu ) dΩ
tfeop = TransientFEOperator(res, (jac, jac_t), U, V)

ls = LUSolver()
nls = NLSolver(ls, show_trace=true, method=:newton, iterations=10)
odeslvr = ThetaMethod(nls, dt, 0.5)

sol = solve(odeslvr, tfeop, t0, tF, uh0)
for (tn, un) in enumerate(sol)
    # ...
end
```

# Low-level implementation
We now briefly describe the low-level implementation of the ODE module in `Gridap`.

## ODE operators
The `ODEOperator` type represents an ODE according to the description above. It implements the `NonlinearOperator` interface, which enables the computation of residuals and jacobians.

The algebraic equivalent of `TransientFEOperator` is an `ODEOpFromTFEOp`, which is a subtype of `ODEOperator`. Conceptually, `ODEOpFromTFEOp` can be thought of as an assembled `TransientFEOperator`, i.e. it deals with vectors of degrees of freedom. This operator comes with a cache (`ODEOpFromTFEOpCache`) that stores the transient space, its evaluation at the current time step, a cache for the `TransientFEOperator` itself (if any), and the constant forms (if any).

> For now `TransientFEOperator` does not implement the `FEOperator` interface, i.e. it is not possible to evaluate residuals and jacobians directly on it. Rather, they are meant to be evaluated on the `ODEOpFromFEOp`. This is to cut down on the number of conversions between a `TransientCellField` and its vectors of degrees of freedom (one per time derivative). Indeed, when linear forms are constant, no conversion is needed as the jacobian matrix will be stored.

## ODE solvers
An ODE solver has to implement the following interface.
* `allocate_odecache(odeslvr, odeop, t0, us0)`. This function allocates a cache that can be reused across the three functions `ode_start`, `ode_march!`, and `ode_finish!`. In particular, it is necessary to call `allocate_odeopcache` within this function, so as to instantiate the `ODEOpFromTFEOpCache` and be able to update the Dirichlet boundary conditions in the subsequent functions.
* `ode_start(odeslvr, odeop, t0, us0, odecache)`. This function creates the state vectors from the initial conditions. By default, this is the identity.
* `ode_march!(stateF, odeslvr, odeop, t0, state0, odecache)`. This is the update map that evolves the state vectors.
* `ode_finish!(uF, odeslvr, odeop, t0, tF, stateF, odecache)`. This function converts the state vectors into the evaluation of the solution at the current time step. By default, this copies the first state vector into `uF`.

## Stage operator
A `StageOperator` represents the linear or nonlinear operator that a numerical scheme relies on to evolve the state vector. It is essentially a special kind of `NonlinearOperator` but it overwrites the behaviour of nonlinear and linear solvers to take advantage of the matrix of a linear system being constant. The following subtypes of `StageOperator` are the building blocks of all numerical schemes.
* `LinearStageOperator` represents the system ``\boldsymbol{J} \boldsymbol{x} + \boldsymbol{r} = \boldsymbol{0}``, and can build ``\boldsymbol{J}`` and ``\boldsymbol{r}`` by evaluating the residual at a given point.
* `NonlinearStageOperator` represents ``\boldsymbol{r}(\boldsymbol{t}, \boldsymbol{\ell}_{0}(\boldsymbol{x}), \ldots, \boldsymbol{\ell}_{N}(\boldsymbol{x})) = \boldsymbol{0}``, where it is assumed that all the ``\boldsymbol{\ell}_{k}(\boldsymbol{x})`` are linear in ``\boldsymbol{x}``.

## ODE solution
This type is a simple wrapper around an `ODEOperator`, an `ODESolver`, and initial conditions that can be iterated on to evolve the ODE.

# Numerical schemes formulation and implementation
We conclude this note by describing some numerical schemes and their implementation in `Gridap`.

Suppose that the scheme has been evolved up to time ``t_{n}`` already and that the state vectors ``\{\boldsymbol{s}\}_{n}`` are known. We are willing to evolve the ODE up to time ``t_{n+1} > t_{n}``, i.e. compute the state vectors ``\{\boldsymbol{s}\}_{n+1}``. Generally speaking, a numerical scheme constructs an approximation of the map ``\{\boldsymbol{s}\}_{n} \to \{\boldsymbol{s}\}_{n+1}`` by solving one or more relationships of the type
```math
\boldsymbol{r}(t_{i}, \Delta_{i}^{0}(\{\boldsymbol{s}\}_{n}, \{\boldsymbol{s}\}_{n+1}), \ldots, \Delta_{i}^{n}(\{\boldsymbol{s}\}_{n}, \{\boldsymbol{s}\}_{n+1})) = \boldsymbol{0},
```
where ``t_{i}`` is an intermediate time and ``\Delta_{i}^{k}`` are discrete operators that approximates the ``k``-th order time derivative of ``\boldsymbol{u}`` at time ``t_{i}``.

We now describe the numerical schemes implemented in `Gridap` using this framework. It is usually convenient to perform a change of variables so that the unknown ``\boldsymbol{x}`` has the dimension of the highest-order time derivative of ``\boldsymbol{u}``, i.e. ``[\boldsymbol{x}] = [t]^{-n} [\boldsymbol{u}]`` (where ``[\bullet]`` stands for "the dimension of ``\bullet``"). We always perform such a change of variable in practice.

We also briefly characterise these schemes in terms of their order and linear stability.

## ``\theta``-method
This scheme is used to solve first-order ODEs and relies on the simple state vector ``\{\boldsymbol{s}(t)\} = \{\boldsymbol{u}(t)\}``. This means that the starting and finishing procedures are simply the identity.

The ``\theta``-method relies on the following approximation
```math
\boldsymbol{u}(t_{n+1}) = \boldsymbol{u}(t_{n}) + \int_{t_{n}}^{t_{n+1}} \partial_{t} \boldsymbol{u}(t) \ \mathrm{d} t \approx \boldsymbol{u}(t_{n}) + h_{n} \partial_{t} \boldsymbol{u}(t_{n + \theta}),
```
where we have introduced the intermediate time ``t_{n + \theta} \doteq (1 - \theta) t_{n} + \theta t_{n+1}``. By replacing ``\boldsymbol{u}(t_{n})`` and ``\boldsymbol{u}(t_{n+1})`` by their discrete equivalents, we have ``\partial_{t} \boldsymbol{u}(t_{n + \theta}) \approx \frac{1}{h} (\boldsymbol{u}_{n+1} - \boldsymbol{u}_{n})``. This quantity is found by enforcing that the residual is zero at ``t_{n + \theta}``. In that sense, the ``\theta``-method can be framed as a collocation method at ``t_{n + \theta}``. For that purpose, we use the same quadrature rule as above to approximate ``\boldsymbol{u}(t_{n + \theta})``, i.e. ``\boldsymbol{u}(t_{n + \theta}) \approx \boldsymbol{u}_{n} + \theta h_{n} \partial_{t} \boldsymbol{u}(t_{n + \theta})``. Using the notations of the framework above, we have defined
```math
\begin{align*}
t_{1} &= (1 - \theta) t_{n} + \theta t_{n+1}, \\
\Delta_{1}^{0} &= (1 - \theta) \boldsymbol{u}_{n} + \theta \boldsymbol{u}_{n+1}, \\
\Delta_{1}^{1} &= \frac{1}{h} (\boldsymbol{u}_{n+1} - \boldsymbol{u}_{n}).
\end{align*}
```

To summarize and to be more concrete, let ``\boldsymbol{x} = \frac{1}{h} (\boldsymbol{u}_{n+1} - \boldsymbol{u}_{n})``. The ``\theta``-method solves the following stage operator
```math
\boldsymbol{r}(t_{n} + \theta h_{n}, \boldsymbol{u}_{n} + \theta h_{n} \boldsymbol{x}, \boldsymbol{x}) = \boldsymbol{0},
```
and sets ``\boldsymbol{u}_{n+1} = \boldsymbol{u}_{n} + h_{n} \boldsymbol{x}``. The output state is simply ``\{\boldsymbol{s}\}_{n+1} = \{\boldsymbol{u}_{n+1}\}``.

##### Analysis
Since this scheme uses ``\boldsymbol{u}(t)`` as its only state vector, the amplification matrix has dimension one, and its coefficient is the stabilisation function, given by
```math
\rho(z) = \frac{1 + (1 - \theta) z}{1 - \theta z}.
```
We plug the Taylor expansion of ``\boldsymbol{u}_{n+1}`` around ``\boldsymbol{u}_{n}`` in ``\boldsymbol{u}_{n+1} = \rho(z) \boldsymbol{u}_{n}`` and obtain the exactness condition ``\rho(z) - \exp(z) = 0``. We then seek to match as many coefficients in the Taylor expansion of both sides to obtain order conditions. We readily obtain the following expansion
```math
\rho(z) - \exp(z) = \sum_{k \geq 0} \left[\theta^{k} - \frac{1}{(k+1)!}\right] z^{k+1}.
```
The order conditions are as follows.
* **Order 0 and 1**. The first two coefficients are always zero, so the method has at least order one.
* **Order 2**. The third coefficient is ``\theta - \frac{1}{2}``, and it is zero when ``\theta = \frac{1}{2}``. This value of ``\theta`` corresponds to a second-order scheme. The next coefficient is ``\theta^{2} - \frac{1}{6}``, so this method cannot reach order three.

By looking at the behaviour of the stability function at infinity, we find that the scheme is ``L``-stable only when ``\theta = 1``. We determine whether the scheme is ``A``-stable or not by looking at stability region. We distinguish three cases based on the value of ``\theta``.
* ``\theta < \frac{1}{2}``. The stability region is the circle of radius ``\frac{1}{1 - 2 \theta}`` centered at ``\left(\frac{-1}{1 - 2 \theta}, 0\right)``. In particular, it is not ``A``-stable. The special case ``\theta = 0`` is known as the Forward Euler scheme, which is the only explicit scheme of the ``\theta``-method family.
* ``\theta = \frac{1}{2}``. The stability region is the whole left complex plane, so the scheme is ``A``-stable. This case is known as the implicit midpoint scheme. 
* ``\theta > \frac{1}{2}``. The stability region is the whole complex plane except the circle of radius ``\frac{1}{2 \theta - 1}`` centered at ``\left(\frac{1}{2 \theta - 1}, 0\right)``. In particular, the scheme is ``A``-stable. The special case ``\theta = 1`` is known as the Backward Euler scheme. 

## Generalised-``\alpha`` scheme for first-order ODEs
This scheme relies on the state vector ``\{\boldsymbol{s}(t)\} = \{\boldsymbol{u}(t), \partial_{t} \boldsymbol{u}(t)\}``. In particular, it needs a nontrivial starting procedure that evaluates ``\partial_{t} \boldsymbol{u}(t_{0})`` by enforcing a zero residual at ``t_{0}``. The finaliser can still return the first vector of the state vectors. For convenience, let ``\partial_{t} \boldsymbol{u}_{n}`` denote the approximation ``\partial_{t} \boldsymbol{u}(t_{n})``.

> Alternatively, the initial velocity can be provided manually: when calling `solve(odeslvr, tfeop, t0, tF, uhs0)`, set `uhs0 = (u0, v0, a0)` instead of `uhs0 = (u0, v0)`. This is useful when enforcing a zero initial residual would lead to a singular system.

This method extends the ``\theta``-method by considering the two-point quadrature rule
```math
\boldsymbol{u}(t_{n+1}) = \boldsymbol{u}_{n} + \int_{t_{n}}^{t_{n+1}} \partial_{t} \boldsymbol{u}(t) \ \mathrm{d} t \approx \boldsymbol{u}_{n} + h_{n} [(1 - \gamma) \partial_{t} \boldsymbol{u}(t_{n}) + \gamma \partial_{t} \boldsymbol{u}(t_{n+1})],
```
where ``0 \leq \gamma \leq 1`` is a free parameter. The question is now how to estimate ``\partial_{t} \boldsymbol{u}(t_{n+1})``. This is achieved by enforcing a zero residual at ``t_{n + \alpha_{F}} \doteq (1 - \alpha_{F}) t_{n} + \alpha_{F} t_{n+1}``, where ``0 \leq \alpha_{F} \leq 1`` is another free parameter. The value of ``\boldsymbol{u}`` at that time, ``\boldsymbol{u}_{n + \alpha_{F}}``, is obtained by the same linear combination of ``\boldsymbol{u}`` at ``t_{n}`` and ``t_{n+1}``. Regarding ``\partial_{t} \boldsymbol{u}``, it is taken as a linear combination weighted by another free parameter ``0 < \alpha_{M} \leq 1`` of the time derivative at times ``t_{n}`` and ``t_{n+1}``. Note that ``\alpha_{M}`` cannot be zero. Altogether, we have defined the discrete operators
```math
\begin{align*}
t_{1} &= (1 - \alpha_{F}) t_{n} +  \alpha_{F} t_{n+1}, \\
\Delta_{1}^{0} &= (1 - \alpha_{F}) \boldsymbol{u}_{n} + \alpha_{F} \boldsymbol{u}_{n+1}, \\
\Delta_{1}^{1} &= (1 - \alpha_{M}) \partial_{t} \boldsymbol{u}_{n} + \alpha_{M} \partial_{t} \boldsymbol{u}_{n+1}.
\end{align*}
```

In more concrete terms, we solve the following system:
```math
\begin{align*}
\boldsymbol{0} &= \boldsymbol{r}(t_{n + \alpha_{F}}, \boldsymbol{u}_{n + \alpha_{F}}, \partial_{t} \boldsymbol{u}_{n + \alpha_{M}}), \\
t_{n + \alpha_{F}} &= (1 - \alpha_{F}) t_{n} + \alpha_{F} t_{n+1}, \\
\boldsymbol{u}_{n + \alpha_{F}} &= (1 - \alpha_{F}) \boldsymbol{u}_{n} + \alpha_{F} \boldsymbol{u}_{n+1}, \\
\partial_{t} \boldsymbol{u}_{n + \alpha_{M}} &= (1 - \alpha_{M}) \partial_{t} \boldsymbol{u}_{n} + \alpha_{M} \partial_{t} \boldsymbol{u}_{n+1}, \\
\boldsymbol{u}_{n+1} &= \boldsymbol{u}_{n} + h_{n} [(1 - \gamma) \partial_{t} \boldsymbol{u}_{n} + \gamma \boldsymbol{x}], \\
\partial_{t} \boldsymbol{u}_{n+1} &= \boldsymbol{x}.
\end{align*}
```
The state vector is updated to ``\{\boldsymbol{s}\}_{n+1} = \{\boldsymbol{u}_{n+1}, \partial_{t} \boldsymbol{u}_{n+1}\}``.

##### Analysis
The amplification matrix for the state vector is
```math
\boldsymbol{A}(z) = \frac{1}{\alpha_{M} - \alpha_{F} \gamma z} \begin{bmatrix}\alpha_{M} + (1 - \alpha_{F}) \gamma z & \alpha_{M} - \gamma \\ z & \alpha_{M} - 1 + \alpha_{F} (1 - \gamma) z\end{bmatrix}.
```
It is then immediate to see that ``\boldsymbol{u}_{n+1} = \mathrm{tr}(\boldsymbol{A}) \boldsymbol{u}_{n} - \det(\boldsymbol{A}) \boldsymbol{u}_{n-1}``. This time, plugging the Taylor expansion of ``\boldsymbol{u}_{n+1}`` and ``\boldsymbol{u}_{n-1}`` around ``\boldsymbol{u}_{n}`` in this expression, the exactness condition is ``\mathrm{tr}(\boldsymbol{A}(z)) - \det(\boldsymbol{A}(z)) \exp(-z) - \exp(z) = 0``. To simplify the analysis, we write the trace and determinant of ``\boldsymbol{A}`` as follows
```math
\mathrm{tr}(\boldsymbol{A}(z)) = a + \frac{b}{1 - c z}, \qquad \det(\boldsymbol{A}(z)) = d + \frac{e}{1 - c z},
```
where
```math
\begin{align*}
a &= 2 - \frac{1}{\alpha_{F}} - \frac{1}{\gamma}, \\
b &= \frac{1}{\alpha_{F}} + \frac{1}{\gamma} - \frac{1}{\alpha_{M}}, \\
c &= \frac{\alpha_{F} \gamma}{\alpha_{M}}, \\
d &= \frac{(1 - \alpha_{F}) (1 - \gamma)}{\alpha_{F} \gamma}, \\
e &= \frac{\alpha_{M} (\alpha_{F} + \gamma - 1) - \alpha_{F} \gamma}{\alpha_{F} \alpha_{M} \gamma}.
\end{align*}
```
Next, we obtain the Taylor expansion of the exactness condition and find
```math
(a + b - d - e - 1) + \sum_{k \geq 1} \left(b c^{k} - \frac{1}{k!} - \frac{(-1)^{k}}{k!} d - \sum_{0 \leq l \leq k} e c^{(k - l)}\frac{(-1)^{l}}{l!}\right) z^{k} = 0.
```
The order conditions are as follows.
* **Order 0 and 1**. The first two coefficients are always zero, so the method is at least of order ``2``.
* **Order 2**. The third coefficient has a zero at ``\gamma = \frac{1}{2} + \alpha_{M} - \alpha_{F}``.
* **Order 3**. The fourth coefficient has a zero at ``\alpha_{M} = \frac{1 + 6 \alpha_{F} - 12 \alpha_{F}^{2}}{6(1 - 2 \alpha_{F})}`` (provided that ``\alpha_{F} \neq \frac{1}{2}``). In that case we simplify ``\gamma`` into ``\gamma = \frac{2 - 3 \alpha_{F}}{3(1 - 2 \alpha_{F})}``.
* **Order 4**. The fifth coefficient has zeros at ``\alpha_{F} = \frac{3 \pm \sqrt{3}}{6}`` and poles at ``\alpha_{F} = \frac{3 \pm \sqrt{21}}{12}``. The corresponding values of ``\alpha_{M}`` and ``\gamma`` are ``\alpha_{M} = \frac{1}{2}``, ``\gamma = \frac{3 \mp \sqrt{3}}{6}``.

We finally study the stability in the extreme cases ``|z| \to 0`` and ``|z| \to +\infty``. We want the spectral radius of the amplification matrix to be smaller than one so that perturbations are damped away.
* When ``|z| \to 0``, we have ``\rho(\boldsymbol{A}(z)) \to \max\{1, \left|1 - \frac{1}{\alpha_{M}}\right|\}``.
* When ``|z| \to +\infty``, we have ``\rho(\boldsymbol{A}(z)) \to \max\{\left|1 - \frac{1}{\alpha_{F}}\right|, \left|1 - \frac{1}{\gamma}\right|\}``.

We thus require ``\alpha_{M} \geq \frac{1}{2}``, ``\alpha_{F} \geq \frac{1}{2}`` and ``\gamma \geq \frac{1}{2}`` to ensure stability. In particular when the scheme has order ``3``, the stability conditions become ``\alpha_{M} \geq \alpha_{F} \geq \frac{1}{2}``. We verify that the scheme is unstable whenever it has an order greater than ``3``. We notice that ``L``-stability is only achieved when ``\alpha_{F} = 1`` and ``\gamma = 1``. The corresponding value of ``\alpha_{M}`` for a third-order scheme is ``\alpha_{M} = \frac{3}{2}``.

This scheme was originally devised to control the damping of high frequencies. One parameterisation consists in prescribing the eigenvalues at ``|z| \to +\infty``, and this leads to
```math
\alpha_{F} = \gamma = \frac{1}{1 + \rho_{\infty}}, \qquad \alpha_{M} = \frac{3 - \rho_{\infty}}{2 (1 + \rho_{\infty})},
```
where ``\rho_{\infty}`` is the spectral radius at infinity. Setting ``\rho_{\infty}`` cuts all the highest frequencies in one step, whereas taking ``\rho_{\infty} = 1`` preserves high frequencies.

## Runge-Kutta
Runge-Kutta methods are multi-stage, i.e. they build estimates of ``\boldsymbol{u}`` at intermediate times between ``t_{n}`` and ``t_{n+1}``. They can be written as follows
```math
\begin{align*}
\boldsymbol{0} &= \boldsymbol{r}(t_{n} + c_{i} h_{n} , \boldsymbol{u}_{n} + \sum_{1 \leq j \leq s} a_{ij} h_{n} \boldsymbol{x}_{j}, \boldsymbol{x}_{i}), & 1 \leq i \leq p \\
\boldsymbol{u}_{n+1} &= \boldsymbol{u}_{n} + \sum_{1 \leq i \leq p} b_{i} h_{n} \boldsymbol{x}_{i},
\end{align*}
```
where ``p`` is the number of stages, ``\boldsymbol{A} = (a_{ij})_{1 \leq i, j \leq p}`` is a matrix of free parameters, ``\boldsymbol{b} = (b_{i})_{1 \leq i \leq p}`` and ``\boldsymbol{c} = (c_{i})_{1 \leq i \leq p}`` are two vectors of free parameters. The stage unknowns ``(\boldsymbol{x}_{i})_{1 \leq i \leq p}`` are involved in a coupled system of equations. This system can take a simpler form when the matrix ``\boldsymbol{A}`` has a particular structure.
* When ``\boldsymbol{A}`` is lower triangular, the equations are decoupled and can thus be solved sequentially. These schemes are called Diagonally-Implicit Runge-Kutta (DIRK). If the diagonal coefficients of the matrix ``\boldsymbol{A}`` are the same, the method is called Singly-Diagonally Implicit (SDIRK).
* If the diagonal coefficients are also zero, the method is explicit. These schemes are called Explicit Runge-Kutta (EXRK).

**Implementation details** It is particularly advantageous to save the factorisation of the matrices of the stage operators for Runge-Kutta methods. This is always possible when the method is explicit and the mass matrix is constant, in which case all the stage matrices are the mass matrix. When the method is diagonally-implicit and the stiffness and mass matrices are constant, the matrices of the stage operators are ``\boldsymbol{M} + a_{ii} h_{n} \boldsymbol{K}``. In particular, if two diagonal coefficients coincide, the corresponding operators will have the same matrix. We implement these reuse strategies by storing them in `CompressedArray`s, and introducing a map `i -> NumericalSetup`.

##### Analysis
The stability function of a Runge-Kutta scheme is
```math
\rho(z) = 1 + z \boldsymbol{b}^{T} (\boldsymbol{I} - z \boldsymbol{A})^{-1} \boldsymbol{1}.
```

The analysis of Runge-Kutta methods is well-established but we only derive order conditions for schemes with one, two, or three stages in the diagonally-implicit case.
* **One stage**. These schemes coincide with the ``\theta``-method presented above.

* **Two stages**. We solve the order conditions given by the differential trees and find the following families of tableaus of orders two and three
```math
\def\arraystretch{1.5}
\begin{array}{c|cc}
\alpha & \alpha & \\
\beta & \beta - \hat{\beta} & \hat{\beta} \\ \hline
& \frac{2 \beta - 1}{2 (\beta - \alpha)} & \frac{1 - 2 \alpha}{2 (\beta - \alpha)}
\end{array}, \qquad
\begin{array}{c|cc}
\frac{1}{2} - \frac{\sqrt{3}}{6} \frac{1}{\lambda} & \frac{1}{2} - \frac{\sqrt{3}}{6} \frac{1}{\lambda} & \\
\frac{1}{2} + \frac{\sqrt{3}}{6} \lambda & \frac{\sqrt{3}}{3} \lambda & \frac{1}{2} - \frac{\sqrt{3}}{6} \lambda \\ \hline
& \frac{\lambda^{2}}{\lambda^{2} + 1} & \frac{1}{\lambda^{2} + 1}.
\end{array}
```

* **Three stages**. We only solve the explicit schemes in full generality. We find three families of order three
```math
\def\arraystretch{1.5}
\begin{array}{c|cc}
0 & \\
\alpha & \alpha & \\
\beta & \beta - \frac{\beta (\beta - \alpha)}{\alpha (2 - 3\alpha)} & \frac{\beta (\beta - \alpha)}{ \alpha(2 - 3 \alpha)} \\ \hline
& 1 - \frac{3 (\beta + \alpha) - 2}{6 \alpha \beta} & \frac{3 \beta - 2}{6 \alpha (\beta - \alpha)} & \frac{2 - 3 \alpha}{6 \beta (\beta - \alpha)}
\end{array}, \qquad
\begin{array}{c|cc}
0 & \\
\frac{2}{3} & \frac{2}{3} & \\
\frac{2}{3} & \frac{2}{3} - \frac{1}{4 \alpha} & \frac{1}{4 \alpha} \\ \hline
& \frac{1}{4} & \frac{3}{4} - \alpha & \alpha
\end{array}, \qquad
\begin{array}{c|cc}
0 & \\
\frac{2}{3} & \frac{2}{3} & \\
0 & -\frac{1}{4 \alpha} & \frac{1}{4 \alpha} \\ \hline
& \frac{1}{4} - \alpha & \frac{3}{4} & \alpha
\end{array}.
```

## Implicit-Explicit Runge-Kutta
When the residual has an implicit-explicit decomposition, usually because we can identify a stiff part that we want to solve implicitly and a nonstiff part that we want to solve explicitly, the Runge-Kutta method reads as follows
```math
\begin{align*}
\boldsymbol{0} &= \boldsymbol{r}(t_{n} + c_{i} h_{n}, \boldsymbol{u}_{n} + \sum_{1 \leq j \leq i-1} (a_{i, j} h_{n} \boldsymbol{x}_{j} + \hat{a}_{i, j} h_{n} \hat{\boldsymbol{x}}_{j}) + a_{i, i} h_{n} \boldsymbol{x}_{i}, \boldsymbol{x}_{i}), \\
\boldsymbol{0} &= \hat{\boldsymbol{r}}(t_{n} + c_{i} h_{n}, \boldsymbol{u}_{n} + \sum_{1 \leq j \leq i-1} (a_{i, j} h_{n} \boldsymbol{x}_{j} + \hat{a}_{i, j} h_{n} \hat{\boldsymbol{x}}_{j}) + a_{i, i} h_{n} \boldsymbol{x}_{i}, \hat{\boldsymbol{x}}_{i}), & 1 \leq i \leq p \\
\boldsymbol{u}_{n+1} &= \boldsymbol{u}_{n} + \sum_{1 \leq i \leq p} (b_{i} h_{n} \boldsymbol{x}_{i} + \hat{b}_{i} h_{n} \hat{\boldsymbol{x}}_{i}).
\end{align*}
```
In these expressions, quantities that wear a hat are the explicit counterparts of the implicit quantity with the same name. The implicit and explicit stages are alternated, i.e. the implicit and explicit stage unknowns ``\boldsymbol{x}_{i}`` and ``\hat{\boldsymbol{x}}_{i}`` are solved alternatively. As seen above, we require that the nodes ``c_{i}`` of the implicit and explicit tableaus coincide. This implies that the first step for the implicit part is actually explicit.

**Implementation details**
Many methods can be created by padding a DIRK tableau with zeros to give it an additional step. In this case, the first stage for the implicit part does not need to be solved, as all linear combinations give it a zero weight. As an example, an ``L``-stable, ``2``-stage, second-order SDIRK IMEX scheme is given by
```math
\def\arraystretch{1.5}
\begin{array}{c|ccc}
0 & 0 & & \\
\frac{2 - \sqrt{2}}{2} & 0 & \frac{2 - \sqrt{2}}{2} & \\
1 & 0 & \frac{\sqrt{2}}{2} & \frac{2 - \sqrt{2}}{2} \\ \hline
 & 0 & \frac{\sqrt{2}}{2} & \frac{2 - \sqrt{2}}{2}
\end{array}, \qquad
\begin{array}{c|ccc}
0 & & & \\
\frac{2 - \sqrt{2}}{2} & \frac{2 - \sqrt{2}}{2} & & \\
1 & -\frac{\sqrt{2}}{2} & 1 + \frac{\sqrt{2}}{2} & \\ \hline
 & -\frac{\sqrt{2}}{2} & 1 + \frac{\sqrt{2}}{2} &
\end{array}.
```
We note that the first column of the matrix and the first weight are all zero, so the first stage for the implicit part does not need to be solved.

## Generalised-``\alpha`` scheme for second-order ODEs
This scheme relies on the state vector ``\{\boldsymbol{s}(t)\} = \{\boldsymbol{u}(t), \partial_{t} \boldsymbol{u}(t), \partial_{tt} \boldsymbol{u}(t)\}``. It needs a nontrivial starting procedure that evaluates ``\partial_{tt} \boldsymbol{u}(t_{0})`` by enforcing a zero residual at ``t_{0}``. The finaliser can still return the first vector of the state vectors. For convenience, let ``\partial_{tt} \boldsymbol{u}_{n}`` denote the approximation ``\partial_{tt} \boldsymbol{u}(t_{n})``.

> The initial acceleration can alternatively be provided manually: when calling `solve(odeslvr, tfeop, t0, tF, uhs0)`, set `uhs0 = (u0, v0, a0)` instead of `uhs0 = (u0, v0)`. This is useful when enforcing a zero initial residual would lead to a singular system.

This method is built out of the following update rule
```math
\begin{align*}
\boldsymbol{0} &= \boldsymbol{r}(t_{n + 1 - \alpha_{F}}, \boldsymbol{u}_{n + 1 - \alpha_{F}}, \partial_{t} \boldsymbol{u}_{n + 1 - \alpha_{F}}, \partial_{tt} \boldsymbol{u}_{n + 1 - \alpha_{M}}), \\
t_{n + 1 - \alpha_{F}} &= \alpha_{F} t_{n} + (1 - \alpha_{F}) t_{n+1}, \\
\boldsymbol{u}_{n + 1 - \alpha_{F}} &= \alpha_{F} \boldsymbol{u}_{n} + (1 - \alpha_{F}) \boldsymbol{u}_{n+1}, \\
\partial_{t} \boldsymbol{u}_{n + 1 - \alpha_{F}} &= \alpha_{F} \partial_{t} \boldsymbol{u}_{n} + (1 - \alpha_{F}) \partial_{t} \boldsymbol{u}_{n+1}, \\
\partial_{tt} \boldsymbol{u}_{n + 1 - \alpha_{M}} &= \alpha_{M} \partial_{tt} \boldsymbol{u}_{n} + (1 - \alpha_{M}) \partial_{tt} \boldsymbol{u}_{n+1}, \\
\boldsymbol{u}_{n+1} &= \boldsymbol{u}_{n} + h_{n} \partial_{t} \boldsymbol{u}_{n} + \frac{1}{2} h_{n}^{2} [(1 - 2 \beta) \partial_{tt} \boldsymbol{u}_{n} + 2 \beta \boldsymbol{x}] \\
\partial_{t} \boldsymbol{u}_{n+1} &= \partial_{t} \boldsymbol{u}_{n} + h_{n} [(1 - \gamma) \partial_{tt} \boldsymbol{u}_{n} + \gamma \boldsymbol{x}], \\
\partial_{tt} \boldsymbol{u}_{n+1} &= \boldsymbol{x}
\end{align*}
```
The state vector is then updated to ``\{\boldsymbol{s}\}_{n+1} = \{\boldsymbol{u}_{n+1}, \partial_{t} \boldsymbol{u}_{n+1}, \partial_{tt} \boldsymbol{u}_{n+1}\}``.

##### Analysis
The amplification matrix for the state vector is
```math
\boldsymbol{A}(z) = \frac{1}{\overline{\alpha_{M}} + \overline{\alpha_{F}} \beta z^{2}} \begin{bmatrix}
\overline{\alpha_{M}} - \alpha_{F} \beta z^{2} & \overline{\alpha_{M}} & \overline{\beta} \overline{\alpha_{M}} - \beta \alpha_{M} \\
-\gamma z^{2} & \overline{\alpha_{M}} + \overline{\alpha_{F}} (\beta - \gamma) z^{2} & \overline{\alpha_{M}} \ \overline{\gamma} - \alpha_{M} \gamma + \overline{\alpha_{F}} [\overline{\gamma} \beta - \overline{\beta} \gamma] z^{2} \\
-z^{2} & -\overline{\alpha_{F}} z^{2} & -\alpha_{M} - \overline{\alpha_{F}} \overline{\beta} z^{2}
\end{bmatrix},
```
where ``\overline{\alpha_{M}} = 1 - \alpha_{M}``, ``\overline{\alpha_{F}} = 1 - \alpha_{F}``, ``\overline{\gamma} = 1 - \gamma`` and ``\overline{\beta} = \frac{1}{2}(1 - 2 \beta)``. Here again, we immediately see that ``\boldsymbol{u}_{n+1}`` satisfies the recurrence
```math
\boldsymbol{u}_{n+1} = \mathrm{tr}(\boldsymbol{A}(z)) \boldsymbol{u}_{n} - \frac{1}{2} (\mathrm{tr}(\boldsymbol{A}(z))^{2} - \mathrm{tr}(\boldsymbol{A}(z)^{2})) \boldsymbol{u}_{n-1} + \det(\boldsymbol{A}(z)) \boldsymbol{u}_{n-2}.
```
By plugging the Taylor expansion of ``\boldsymbol{u}`` at times ``t_{n+1}``, ``t_{n-1}`` and ``t_{n-2}``, we obtain the exactness condition
```math
\cos(z) = \mathrm{tr}(\boldsymbol{A}(z)) - \frac{1}{2} (\mathrm{tr}(\boldsymbol{A}(z))^{2} - \mathrm{tr}(\boldsymbol{A}(z)^{2})) \cos(z) + \det(\boldsymbol{A}(z)) \cos(2z).
```
These conditions are hard to examine analytically, but one can verify that this scheme is at least of order ``1``. Second-order is achieved by setting ``\gamma = \frac{1}{2} - \alpha_{M} + \alpha_{F}``.

It is easier to consider the limit cases ``|z| \to 0`` and ``|z| \to +\infty`` and look at the eigenvalues of the amplification matrix.
* When ``|z| \to 0``, we find ``\rho(\boldsymbol{A}(z)) = \max\left\{1, \left|\frac{\alpha_{M}}{1 - \alpha_{M}}\right|\right\}``.
* When ``|z| \to +\infty``, we find ``\rho(\boldsymbol{A}(z)) = \max\left\{\left|\frac{\alpha_{F}}{1 - \alpha_{F}}\right|, \left|\frac{4 \beta - (1 + 2 \gamma) \pm \sqrt{(1 + 2 \gamma)^{2} - 16 \beta}}{4 \beta}\right|\right\}``.

For all these eigenvalues to have a modulus smaller than one, we need ``\alpha_{M} \leq \frac{1}{2}``, ``\alpha_{F} \leq \frac{1}{2}``, ``\gamma \geq \frac{1}{2}``, i.e. ``\alpha_{F} \geq \alpha_{M}`` and ``\beta \geq \frac{1}{2} \gamma``. Since dissipation of high-frequency is maximised when the eigenvalues are real at infinity, we also impose ``\beta = \frac{1}{16} (1 + 2 \gamma)^{2}``, i.e. ``\beta = \frac{1}{4} (1 - \alpha_{M} + \alpha_{F})^{2}``.

This method was also designed to damp high-frequency perturbations so it is common practice to parameter this scheme in terms of its spectral radius.
* The Hilbert-Huges-Taylor-``\alpha`` (HHT-``\alpha``) method is obtained by setting ``\alpha_{M} = 0``, ``\alpha_{F} = \frac{1 - \rho_{\infty}}{1 + \rho_{\infty}}``.
* The Wood-Bossak-Zienkiewicz-``\alpha`` (WBZ-``\alpha``) method is recovered by setting ``\alpha_{F} = 0`` and ``\alpha_{M} = \frac{\rho_{\infty} - 1}{\rho_{\infty} + 1}``.
* The standard generalised-``\alpha`` method is obtained by setting ``\alpha_{M} = \frac{2 \rho_{\infty - 1}}{\rho_{\infty} + 1}``, ``\alpha_{F} = \frac{\rho_{\infty}}{\rho_{\infty} + 1}``.
* The Newmark method corresponds to ``\alpha_{F} = \alpha_{M} = 0``. In this case, the values of ``\beta`` and ``\gamma`` are usually chosen as ``\beta = 0``, ``\gamma = \frac{1}{2}`` (explicit central difference scheme), or ``\beta = \frac{1}{4}`` and ``\gamma = \frac{1}{2}`` (midpoint rule).

# Reference

```@autodocs
Modules = [ODEs,]
```
