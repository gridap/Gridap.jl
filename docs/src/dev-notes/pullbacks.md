
# Pullbacks

Consider a reference polytope ``\hat{K}``, mapped to the physical space by a **geometrical map** ``G``, i.e. ``K = G(\hat{K})``. Consider also a linear function space on the reference polytope ``\hat{V}``, and a set of unisolvent degrees of freedom represented by moments in the dual space ``\hat{V}^*``.

Throughout this document, we will use the following notation:

- ``\varphi \in V`` is a **physical field** ``\varphi : K \rightarrow \mathbb{R}^k``. A basis of ``V`` is denoted by ``\Phi = \{\varphi\}``.
- ``\hat{\varphi} \in \hat{V}`` is a **reference field** ``\hat{\varphi} : \hat{K} \rightarrow \mathbb{R}^k``. A basis of ``\hat{V}`` is denoted by ``\hat{\Phi} = \{\hat{\varphi}\}``.
- ``\sigma \in V^*`` is a **physical moment** ``\sigma : V \rightarrow \mathbb{R}``. A basis of ``V^*`` is denoted by ``\Sigma = \{\sigma\}``.
- ``\hat{\sigma} \in \hat{V}^*`` is a **reference moment** ``\hat{\sigma} : \hat{V} \rightarrow \mathbb{R}``. A basis of ``\hat{V}^*`` is denoted by ``\hat{\Sigma} = \{\hat{\sigma}\}``.

We define a **pushforward** map as ``F^* : \hat{V} \rightarrow V``, mapping reference fields to physical fields. Given a pushforward ``F^*``, we define:

- The **pullback** ``F_* : V^* \rightarrow \hat{V}^*``, mapping physical moments to reference moments. Its action on physical dofs is defined in terms of the pushforward map ``F^*`` as ``\hat{\sigma} = F_*(\sigma) := \sigma \circ F^*``.
- The **inverse pushforward** ``(F^*)^{-1} : V \rightarrow \hat{V}``, mapping physical fields to reference fields.
- The **inverse pullback** ``(F_*)^{-1} : \hat{V}^* \rightarrow V^*``, mapping reference moments to physical moments. Its action on reference dofs is defined in terms of the inverse pushforward map ``(F^*)^{-1}`` as ``\sigma = (F_*)^{-1}(\hat{\sigma}) := \hat{\sigma} \circ (F^*)^{-1}``.

In many occasions, we will have that (as a basis)

```math
\hat{\Sigma} \neq F_*(\Sigma), \quad \text{and} \quad \Phi \neq F^*(\hat{\Phi})
```

To maintain conformity and proper scaling in these cases, we define cell-dependent invertible changes of basis ``P`` and ``M``, such that

```math
\hat{\Sigma} = P F_*(\Sigma), \quad \text{and} \quad \Phi = M F^*(\hat{\Phi})
```

An important result from linear algebra is that ``P = M^*`` (or just the transpose if we represent them as real matrices).
From an implementation point of view, it is more natural to build ``P^{-1}``  and then retrieve all other matrices by transposition/inversion.

## References

- [Kirby 2017, A general approach to transforming finite elements.](https://arxiv.org/abs/1706.09017)
- [Aznaran et al. 2021, Transformations for Piola-mapped elements.](https://arxiv.org/abs/2110.13224)
