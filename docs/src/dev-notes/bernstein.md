```@meta
CurrentModule = Gridap.Polynomials
```

# Bernstein bases algorithms

### Barycentric coordinates

A ``D``-dimensional simplex ``T`` is defined by ``N=D+1`` vertices ``\{v_1,
v_2, …, v_N\}=\{v_i\}_{i∈1:N}``. The barycentric coordinates
``λ(\bm{x})=\{λ_j(\bm{x})\}_{1 ≤ j ≤ N}`` are uniquely
defined by:
```math
\bm{x} = ∑_{1 ≤ j ≤ N} λ_j(\bm{x})v_j \quad\text{and}\quad
∑_{1≤ j≤ N} λ_j(\bm{x}) = 1,
```
as long as the simplex is non-degenerate (vertices are not all in one
hyperplane).

Assuming the simplex polytopal (has flat faces), this change of coordinates is
affine, and is implemented using:
```math
λ(\bm{x}) = M\left(\begin{array}{c} 1\\ x_1\\ ⋮\\ x_D \end{array}\right)
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
``∂_{x_i} λ_j(\bm{x}) = M_{j,i+1}``, so
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
λ(\bm{x}) = \Big(1-∑_{1≤ i≤ D} x_i, x_1, x_2, ⋯, x_D\Big)
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
B^{D,K}_α(\bm{x}) = \binom{K}{α} λ(\bm{x})^α\qquad\text{for all }α ∈\mathcal{I}_K^D
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
``{\bm{x}_α = \underset{1≤i≤N}{∑} α_i v_i}`` and the intermediate
coefficients used in the de Casteljau algorithm.

These indices are returned by [`bernstein_terms(K,D)`](@ref bernstein_terms).
When storing such quantities in arrays, ``∙_α`` is stored at index
[`bernstein_term_id(α)`](@ref bernstein_term_id), which is the index of `α`
in `bernstein_terms(sum(α),length(α)-1)`.

We adopt the convention that a quantity indexed by a ``α ∉ ℤ_+^N`` is equal to
zero (to simplify the definition of algorithms where ``α=β-e_i`` appears).

### The de Casteljau algorithms

A polynomial ``p ∈ ℙ^D_K`` in Bernstein form ``p = ∑_{α∈\mathcal{I}^D_K}\, p_α
B_α`` can be evaluated at ``\bm{x}`` using the de Casteljau algorithms
[1, Algo. 2.9] by iteratively computing
```math
\qquad p_β^{(l)} = \underset{1 ≤ i ≤ N}{∑} λ_i\, p_{β+e_i}^{(l-1)} \qquad ∀β ∈ \mathcal{I}^D_{K-l},
```
for ``l=1, 2, …, K`` where ``p_α^{(0)}=p_α``, ``λ=λ(\bm{x})`` and the
result is ``p(\bm{x})=p_𝟎^{(K)}``. This algorithm is implemented (in
place) by [`_de_Casteljau_nD!`](@ref Polynomials._de_Casteljau_nD!).

But Gridap implements the polynomial bases themselves instead of individual
polynomials in a basis. To compute all ``B_α`` at ``\bm{x}``, one can
use the de Casteljau algorithm going "downwards" (from the tip of the pyramid
to the base). The idea is to use the relation
```math
B_α = ∑_{1 ≤ i ≤ N} λ_i B_{α-e_i}\qquad ∀α ∈ ℤ_+^N,\ |α|≥1.
```

Starting from ``b_𝟎^{(0)}=B_𝟎(\bm{x})=1``, compute iteratively
```math
\qquad b_β^{(l)} = \underset{1 ≤ i ≤ N}{∑} λ_i\, b_{β-e_i}^{(l-1)} \qquad ∀β ∈ \mathcal{I}^D_{l},
```
for ``l=1,2, …, K``, where again ``λ=λ(\bm{x})`` and the result is
``B_α(\bm{x})=b_α^{(K)}`` for all ``α`` in ``\mathcal{I}^D_K``. This
algorithm is implemented (in place) by [`_downwards_de_Casteljau_nD!`](@ref).
The implementation is a bit tricky, because the iterations must be done in
reverse order to avoid erasing coefficients needed later, and a lot of summands
disappear (when ``(β-e_i)_i < 0``).

The gradient and hessian of the `BernsteinBasisOnSimplex` are also implemented.
They rely on the following
```math
∂_q B_α(\bm{x}) = K\!∑_{1 ≤ i ≤ N} ∂_qλ_i\, B_{α-e_i}(\bm{x}),\qquad
∂_t ∂_q B_α(\bm{x}) = K\!∑_{1 ≤ i,j ≤ N} ∂_tλ_j\, ∂_qλ_i\, B_{α-e_i-e_j}(\bm{x}).
```
The gradient formula comes from [1, Eq. (2.28)], and the second is derived from
the first using the fact that ``∂_qλ`` is homogeneous. The implementation of
the gradient and hessian compute the ``B_β`` using
`_downwards_de_Casteljau_nD!` up to order ``K-1`` and ``K-2`` respectively, and
then the results are assembled by `_grad_Bα_from_Bαm!` and
`_hess_Bα_from_Bαmm!` respectively. The implementation makes sure to only
access each relevant ``B_β`` once per ``(∇/H)B_α`` computed. Also, on the
reference simplex, the barycentric coordinates derivatives are computed at
compile time using ``∂_qλ_i = δ_{i q}-δ_{i N}``.

## Low level docstrings

```@docs
_de_Casteljau_nD!
_downwards_de_Casteljau_nD!
```

# Bernstein basis generalization for ``P_r^{(-)}Λ^k`` spaces

The TODO basis implements the polynomial bases for the spaces ``ℙ_r^-Λ^k`` and
``ℙ_rΛ^k`` (writing ``P_r^{(-)}Λ^k`` for either one of them) on simplices of
any dimension derived in [2]. In this section, we give the translation of the
formulas in the paper from the differential geometry language to usual
differential calculus and the implemented algorithm.

A ``D``-dimensional simplex ``T`` is defined by ``N=D+1`` vertices ``\{v_0, v_1,
..., v_D\}=\{v_i\}_{i\in 1:D}``. We uniquely identify a ``d``-dimensional face
``f`` of ``T`` by the set of the ``d+1`` increasing indices of its vertices:
```math
f\sim I_f = \{i_0, i_1, ..., i_d\} \qquad\text{such that } 0≤ i_0 < i_1 < ... <i_d≤ D .
```
In particular, ``T\sim \{0:D\}``. We write ``f\subset T`` for any face of
``T``, including ``T`` itself or its vertices.
``T`` has ``\binom{N}{d+1}`` ``d``-dimensional faces, indexed ``\forall 0≤ i_0
< i_1 < ... < i_d ≤ D``. The dimension of a face ``f`` is ``\#I_f``, and we
write "``∀\#J=d``" for all the multi-indices of the ``D``-dimensional faces of
``T``.

Using Einstein's convention of summation on repeated indices, a degree-``k``
dimension-``D`` form ``ω`` can be written in the canonical Cartesian basis as
``ω = ω_I\,\text{d}x^I``, where the basis is
```math
\big\{ \text{d}x^I = \underset{j\in I}{\bigwedge}\text{d}x^{j}
=\text{d}x^{i_1}\wedge ...\wedge \text{d}x^{i_k} \quad\big|\quad I=\{i_1, ...,
i_k\} \text{ for }1≤ i_1 < ... < i_k ≤ D\big\},
```
and ``\{ω_I\}_I\in\mathbb{R}^\binom{D}{k}`` is the vector of coefficients.


#### Geometric decomposition
The main feature of the bases ``P_r^{(-)}Λ^k`` is that each basis polynomial
``ω^{f,α}`` is associated with a face ``f`` of ``T`` (and a domain point
``\boldsymbol{x}_α`` strictly inside ``f``), in the sense that the trace of
``ω^{f,α}`` on another face ``g\subset T`` is zero when ``g`` does not contain
``f``:
```math
f\not\subset g\ \rightarrow\ \text{tr}_g ω^{f,α} = 0, \quad\forall f,g \subseteq T,\ \forall \llbracket α\rrbracket\subseteq I_f, α>0.
```
These basis polynomials ``ω^{f,α}`` are called bubble functions associated to
``f``, the space they span is called ``\mathring{ℙ}_r^{(-)}Λ^k(T,f)``. There
are no bubble functions of degree ``k`` on faces of dimension ``<k``, so the
spaces ``ℙ_r^{(-)}Λ^k(T)`` admit the geometric decomposition:
```math
ℙ_r^{(-)}Λ^k(T) = \underset{f\subset T}{\oplus}\ \mathring{ℙ}_r^{(-)}Λ^k(f)
= \underset{k≤d≤D}{\oplus}\underset{\quad I_f=0≤ i_0 < ... < i_d ≤ D}{\oplus}\ \mathring{ℙ}_r^{(-)}Λ^k(T,f).
```

#### Bubble functions ``\mathring{ℙ}_r^-Λ^k``
```math
\mathring{ℙ}_r^-Λ^k(T,f) = \text{span}\big\{B_α φ^J \quad\big|\quad |α|=r-1,\ \#J=k+1,\ ⟦α⟧∪J=I_f,\ α_i=0 \text{ if } i< \text{min}(J) \big\}\newline
```
where ``B_α`` are the scalar Bernstein polynomials implemented by
[`BernsteinBasisOnSimplex`](@ref), and ``φ^J`` are the Whitney forms
```math
φ^J = \sum_{0≤j≤k} (-1)^{j} λ_j \, \text{d}λ^{J\backslash j} \quad\text{where}\quad
\text{d}λ^{J\backslash j} = \underset{i\in J\backslash \{J_j\} }{\bigwedge}\text{d}λ^{i},
```
``φ^J `` is a ``k``-form of polynomial order ``1``. In a polytopal tetrahedron
``T`` (flat faces), the 1-forms ``\text{d}λ^i`` are homogeneous and coefficients
expression in the canonical basis ``\{\text{d}x^I\}_I`` is
```math
```


#### Bubble functions ``\mathring{ℙ}_rΛ^k``
```math
\mathring{ℙ}_rΛ^k(T,f) = \text{span}\big\{B_α ψ_J^{α,f} \quad\big|\quad |α|=r,\ \#J=k,\ ⟦α⟧∪J=I_f,\ α_i=0 \text{ if } i< \text{min}(I_f \backslash J) \big\}
```

```math
```

## References

[1] [M.J. Lai & L.L. Schumaker, Spline Functions on Triangulations, Chapter 2 - Bernstein–Bézier Methods for Bivariate Polynomials, pp. 18 - 61.](https://doi.org/10.1017/CBO9780511721588.003)

[2] [D.N. Arnold, R.S. Falk & R. Winther, Geometric decompositions and local bases for spaces of finite element differential forms, Computer Methods in Applied Mechanics and Engineering](https://doi.org/10.1016/j.cma.2008.12.017)


