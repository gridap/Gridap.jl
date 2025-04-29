# Bernstein bases algorithms

### Barycentric coordinates

A ``D``-dimensional simplex ``T`` is defined by ``N=D+1`` vertices ``\{v_1,
v_2, …, v_N\}=\{v_i\}_{i∈1:N}``. The barycentric coordinates
``λ(\boldsymbol{x})=\{λ_j(\boldsymbol{x})\}_{1 ≤ j ≤ N}`` are uniquely
defined by:
```math
\boldsymbol{x} = ∑_{1 ≤ j ≤ N} λ_j(\boldsymbol{x})v_j \quad\text{and}\quad
∑_{1≤ j≤ N} λ_j(\boldsymbol{x}) = 1,
```
as long as the simplex is non-degenerate (vertices are not all in one
hyperplane).

Assuming the simplex polytopal (has flat faces), this change of coordinates is
affine, and is implemented using:
```math
λ(\boldsymbol{x}) = M\left(\begin{array}{c} 1\\ x_1\\ ⋮\\ x_D \end{array}\right)
\quad\text{with}\quad M =
\left(\begin{array}{cccc}
1 & 1 & ⋯ & 1 \\
(v_1)_1 & (v_2)_1  & ⋯ & (v_N)_1  \\
⋮  & ⋮   & ⋯ & ⋮   \\
(v_1)_D & (v_2)_D  & ⋯ & (v_N)_D  \\
\end{array}\right)^{-1}
```
where the inverse exists because ``T`` is non-degenerate [1], cf. functions
`_cart_to_bary` and `_compute_cart_to_bary_matrix`. Additionally, we have
``∂_{x_i} λ_j(\boldsymbol{x}) = M_{j,i+1}``, so
```math
∇ λ_j = M_{2:N, j}.
```
The matrix ``M`` is all we need that depends on ``T`` in order to compute
Bernstein polynomials and their derivatives, it is stored in the field
`cart_to_bary_matrix` of [`BernsteinBasisOnSimplex`](@ref), when ``T`` is not
the reference simplex.

On the reference simplex defined by the vertices `get_vertex_coordinates(SEGMENT / TRI / TET⋯)`:
```math
\begin{aligned}
v_1     & = (0\ 0\ ⋯\ 0), \\
v_2     & = (1\ 0\ ⋯\ 0), \\
⋮  &         \\
v_N     & = (0\ ⋯\ 0\ 1),
\end{aligned}
```
the matrix ``M`` is not stored because
```math
λ(\boldsymbol{x}) = \Big(1-∑_{1≤ i≤ D} x_i, x_1, x_2, ⋯, x_D\Big)
\quad\text{and}\quad
∂_{x_i} λ_j = δ_{i+1,j} - δ_{1j} = M_{j,i+1}.
```

### Bernstein polynomials definition

The univariate [`Bernstein`](@ref) polynomials forming a basis of ``ℙ_K``
are defined by
```math
B^K_{n}(x) = \binom{K}{n} x^n (1-x)^{K-n}\qquad\text{ for } 0≤ n≤ K.
```

The ``D``-multivariate Bernstein polynomials of degree ``K`` relative to a
simplex ``T`` are defined by
```math
B^{D,K}_α(x) = \binom{K}{α} λ^α\qquad\text{for all }α ∈\mathcal{I}_K^D
```
where
- ``\mathcal{I}_K^D = \{\ α∈(\mathbb{Z}_+)^{D+1} \quad|\quad |α|=K\ \}``
- ``|α|=∑_{1≤ i≤ N} α_i``
- ``\binom{K}{α} = \frac{K!}{α_1 !α_2 !… α_N!}``
- ``λ`` are the barycentric coordinates relative to ``T`` (defined above)

The superscript ``D`` and ``K`` in ``B^{D,K}_α(x)`` can be omitted because they
are always determined by ``α`` using ``{D=\#(α)-1}`` and ``K=|α|``. The set
``\{B_α\}_{α∈\mathcal{I}_K^D}`` is a basis of ``ℙ^D_K``, implemented by
[`BernsteinBasisOnSimplex`](@ref).

### Bernstein indices and indexing

Working with Bernstein polynomials requires dealing with several quantities
indexed by some ``α ∈ \mathcal{I}_K^D``, the polynomials themselves but also the
coefficients ``c_α`` of a polynomial in the basis, the domain points
``{\boldsymbol{x}_α = \underset{1≤i≤N}{∑} α_i v_i}`` and the intermediate
coefficients used in the de Casteljau algorithm.

These indices are returned by [`bernstein_terms(K,D)`](@ref bernstein_terms).
When storing such quantities in arrays, ``∙_α`` is stored at index
[`bernstein_term_id(α)`](@ref bernstein_term_id), which is the index of `α`
in `bernstein_terms(sum(α),length(α)-1)`.

We adopt the convention that a quantity indexed by a ``α ∉ ℤ_+^N`` is equal to
zero (to simplify the definition of algorithms where ``α=β-e_i`` appears).

### The de Casteljau algorithms

A polynomial ``p ∈ ℙ^D_K`` in Bernstein form ``p = ∑_{α∈\mathcal{I}^D_K}\, p_α
B_α`` can be evaluated at ``\boldsymbol{x}`` using the de Casteljau algorithms
[1, Algo. 2.9] by iteratively computing
```math
\qquad p_β^{(l)} = \underset{1 ≤ i ≤ N}{∑} λ_i\, p_{β+e_i}^{(l-1)} \qquad ∀β ∈ \mathcal{I}^D_{K-l},
```
for ``l=1, 2, …, K`` where ``p_α^{(0)}=p_α``, ``λ=λ(\boldsymbol{x})`` and the
result is ``p(\boldsymbol{x})=p_𝟎^{(K)}``. This algorithm is implemented (in
place) by [`_de_Casteljau_nD!`](@ref).

But Gridap implements the polynomial bases themselves instead of individual
polynomials in a basis. To compute all ``B_α`` at ``\boldsymbol{x}``, one can
use the de Casteljau algorithm going "downwards" (from the tip of the pyramid
to the base). Starting from ``b_𝟎^{(K)}=1``, compute iteratively
```math
\qquad b_β^{(l-1)} = \underset{1 ≤ i ≤ N}{∑} λ_i\, b_{β-e_i}^{(l)} \qquad ∀β ∈ \mathcal{I}^D_{K-l+1},
```
for ``l=K-1, …, 1``, where again ``λ=λ(\boldsymbol{x})`` and the result is
``B_α(\boldsymbol{x})=b_α^{(0)}`` for all ``α`` in ``\mathcal{I}^D_K``. This
algorithm is implemented (in place) by [`_downwards_de_Casteljau_nD!`](@ref).
The implementation is a bit tricky, because the iterations must be done in
reverse order to avoid erasing coefficients needed later, and a lot of summands
disappear (when ``(β-e_i)_i < 0``).

### TODO Evaluation of Bernstein polynomials derivatives

```math
B_α = ∑_{1 ≤ i ≤ N} λ_i B_{α-e_i} = ∑_{(i,β) ∈ sub(α)} λ_i B_β
```

```math
∂_q B_α = K∑_{|β|=|α|-1} B_β (
    δ_{β}^{α-e_q}
   -δ_{β}^{α-e_N}
) = K∑_{(i,β) ∈ sub(α)}  ( … )
```

```math
∂_t ∂_q B_α = K(K-1)∑_{|β|=|α|-2} B_β (
    δ_{β}^{α-e_t-e_q}
   -δ_{β}^{α-e_t-e_N}
   -δ_{β}^{α-e_N-e_q}
   +δ_{β}^{α-e_N-e_N}
) = K(K-1)∑_{(i,β) ∈ sub(α,2)}  ( … )
```

```math
```

## References

[1] [M.J. Lai & L.L. Schumaker, Spline Functions on Triangulations, Chapter 2 - Bernstein–Bézier Methods for Bivariate Polynomials, pp. 18 - 61.](https://doi.org/10.1017/CBO9780511721588.003)
