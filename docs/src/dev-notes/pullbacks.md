
# FE basis transformations

## Notation

Consider a reference polytope ``\hat{K}``, mapped to the physical space by a **geometrical map** ``F``, i.e. ``K = F(\hat{K})``. Consider also a linear function space on the reference polytope ``\hat{V}``, and a set of unisolvent degrees of freedom represented by moments in the dual space ``\hat{V}^*``.

Throughout this document, we will use the following notation:

- ``\varphi \in V`` is a **physical field** ``\varphi : K \rightarrow \mathbb{R}^k``. A basis of ``V`` is denoted by ``\Phi = \{\varphi\}``.
- ``\hat{\varphi} \in \hat{V}`` is a **reference field** ``\hat{\varphi} : \hat{K} \rightarrow \mathbb{R}^k``. A basis of ``\hat{V}`` is denoted by ``\hat{\Phi} = \{\hat{\varphi}\}``.
- ``\sigma \in V^*`` is a **physical moment** ``\sigma : V \rightarrow \mathbb{R}``. A basis of ``V^*`` is denoted by ``\Sigma = \{\sigma\}``.
- ``\hat{\sigma} \in \hat{V}^*`` is a **reference moment** ``\hat{\sigma} : \hat{V} \rightarrow \mathbb{R}``. A basis of ``\hat{V}^*`` is denoted by ``\hat{\Sigma} = \{\hat{\sigma}\}``.

## Pullbacks and Pushforwards

We define a **pushforward** map as ``F^* : \hat{V} \rightarrow V``, mapping reference fields to physical fields. Given a pushforward ``F^*``, we define:

- The **pullback** ``F_* : V^* \rightarrow \hat{V}^*``, mapping physical moments to reference moments. Its action on physical dofs is defined in terms of the pushforward map ``F^*`` as ``\hat{\sigma} = F_*(\sigma) := \sigma \circ F^*``.
- The **inverse pushforward** ``(F^*)^{-1} : V \rightarrow \hat{V}``, mapping physical fields to reference fields.
- The **inverse pullback** ``(F_*)^{-1} : \hat{V}^* \rightarrow V^*``, mapping reference moments to physical moments. Its action on reference dofs is defined in terms of the inverse pushforward map ``(F^*)^{-1}`` as ``\sigma = (F_*)^{-1}(\hat{\sigma}) := \hat{\sigma} \circ (F^*)^{-1}``.

## Change of basis

In many occasions, we will have that (as a basis)

```math
\hat{\Sigma} \neq F_*(\Sigma), \quad \text{and} \quad \Phi \neq F^*(\hat{\Phi})
```

To maintain conformity and proper scaling in these cases, we define cell-dependent invertible changes of basis ``P`` and ``M``, such that

```math
\hat{\Sigma} = P F_*(\Sigma), \quad \text{and} \quad \Phi = M F^*(\hat{\Phi})
```

An important result from [1, Theorem 3.1] is that ``P = M^T``.

!!! details
    [1, Lemma 2.6]: A key ingredient is that given ``M`` a matrix we have ``\Sigma (M \Phi) = \Sigma (\Phi) M^T`` since
    ```math
      [\Sigma (M \Phi)]_{ij} = \sigma_i (M_{jk} \varphi_k) = \sigma_i \varphi_k M_{jk} = [\Sigma (\Phi) M^T]_{ij}
    ```
    where we have used that moments are linear.

We then have the following diagram:

```math
\hat{V}^* \xleftarrow{P} \hat{V}^* \xleftarrow{F_*} V^* \\
\hat{V} \xrightarrow{F^*} V \xrightarrow{P^T} V
```

!!! details
    The above diagram is well defined, since we have
    ```math
    \hat{\Sigma}(\hat{\Phi}) = P F_* (\Sigma)(F^{-*} (P^{-T} \Phi)) = P \Sigma (F^* (F^{-*} P^{-T} \Phi)) = P \Sigma (P^{-T} \Phi) = P \Sigma (\Phi) P^{-1} = Id \\
    \Sigma(\Phi) = F_*^{-1}(P^{-1}\hat{\Sigma})(P^T F^*(\hat{\Phi})) = P^{-1} \hat{\Sigma} (F^{-*}(P^T F^*(\hat{\Phi}))) = P^{-1} \hat{\Sigma} (P^T \hat{\Phi}) = P^{-1} \hat{\Sigma}(\hat{\Phi}) P = Id
    ```

From an implementation point of view, it is more natural to build ``P^{-1}``  and then retrieve all other matrices by transposition/inversion.

## Interpolation

In each cell ``K`` and for ``C_b^k(K)`` the space of functions defined on ``K`` with at least ``k`` bounded derivatives, we define the interpolation operator ``I_K : C_b^k(K) \rightarrow V`` as

```math
I_K(g) = \Sigma(g) \Phi \quad, \quad \Sigma(g) = P^{-1} \hat{\Sigma}(F^{-*}(g))
```

## Implementation notes

!!! note
    In [2], Covariant and Contravariant Piola maps preserve exactly (without any sign change) the normal and tangential components of a vector field.
    I am quite sure that the discrepancy is coming from the fact that the geometrical information in the reference polytope is globaly oriented.
    For instance, the normals ``n`` and ``\hat{n}`` both have the same orientation, i.e ``n = (||\hat{e}||/||e||) (det J) J^{-T} \hat{n}``. Therefore ``\hat{n}`` is not fully local. See [2, Equation 2.11].
    In our case, we will be including the sign change in the transformation matrices, which will include all cell-and-dof-dependent information.

## References

[1] [Kirby 2017, A general approach to transforming finite elements.](https://arxiv.org/abs/1706.09017)

[2] [Aznaran et al. 2021, Transformations for Piola-mapped elements.](https://arxiv.org/abs/2110.13224)
