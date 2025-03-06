
# Bernstein basis algorithms

## Notations

The univariate [`Bernstein`](@ref) polynomials forming a basis of ``ℙ_K``
are defined by
```math
B^K_{n}(x) = \binom{K}{n} x^n (1-x)^{K-n}\qquad\text{ for } 0≤ n≤ K.
```

The ``D``-multivariate Bernstein polynomials of degree ``K`` are defined by
```math
B^{D,K}_\alpha(x) = \binom{K}{\alpha} λ^\alpha\qquad\text{for } |\alpha|=K
```
where
- ``\alpha`` belongs to ``\{ \alpha\in \}`` where ``N = D+1``
- ``\binom{K}{\alpha} = \frac{K!}{\alpha_1 !\alpha_2 !\dots\alpha_N!}`` and ``|\alpha|=\sum_{1≤ i≤ N} \alpha_i``

The superscript ``D`` and ``K`` can be omitted because they are always
determined by ``\alpha`` using ``D=\#(α)-1`` and ``K=|α|``.
The set ``\{B_\alpha\}_{\#(\alpha)=D+1,\,|\alpha|=K}`` is a basis of
``ℙ^D_K``, implemented by [`BernsteinBasisOnSimplex`](@ref).

## Pascal triangle, pyramid and simplex



## Standard de Casteljau algorithm


## Downward de Casteljau algorithm


## Evaluation of Bernstein polynomials and their derivatives

```math
B_\alpha = \sum_{1≤ i≤ N} λ_i B_{α-e_i} = \sum_{(i,β)\in sub(α)} λ_i B_β
```

```math
\partial_q B_\alpha = K\sum_{|\beta|=|\alpha|-1} B_\beta (
    \delta_{\beta}^{\alpha-e_q}
   -\delta_{\beta}^{\alpha-e_N}
) = K\sum_{(i,β)\in sub(α)}  ( \dots )
```

```math
\partial_t\partial_q B_\alpha = K(K-1)\sum_{|\beta|=|\alpha|-2} B_\beta (
    \delta_{\beta}^{\alpha-e_t-e_q}
   -\delta_{\beta}^{\alpha-e_t-e_N}
   -\delta_{\beta}^{\alpha-e_N-e_q}
   +\delta_{\beta}^{\alpha-e_N-e_N}
) = K(K-1)\sum_{(i,β)\in sub(α,2)}  ( \dots )
```

```math
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


