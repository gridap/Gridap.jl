
# Bernstein bases algorithms

### Barycentric coordinates

A ``D``-dimensional simplex ``T`` is defined by ``N=D+1`` vertices ``\{v_1,
v_2, ..., v_N\}=\{v_i\}_{i\in 1:N}``. The barycentric coordinates
``λ(\boldsymbol{x})=\{λ_j(\boldsymbol{x})\}_{1\leq j\leq N}`` are uniquely
defined by:
```math
\boldsymbol{x} = \sum_{1\leq j\leq N} λ_j(\boldsymbol{x})v_j \quad\text{and}\quad
\sum_{1\leq j\leq N} λ_j(\boldsymbol{x}) = 1,
```
as long as the simplex is non-degenerate (vertices are not all in one
hyperplane).

In polytopal simplices (those with flat faces), this change of coordinates is affine,
and is implemented using:
```math
λ(\boldsymbol{x}) = M\left(\begin{array}{c} 1\\ x_1\\ ...\\ x_D \end{array}\right)
\quad\text{with}\quad M =
\left(\begin{array}{cccc}
1 & 1 & ... & 1 \\
(v_1)_1 & (v_2)_1  & ... & (v_N)_1  \\
\vdots  & \vdots   & ... & \vdots   \\
(v_1)_D & (v_2)_D  & ... & (v_N)_D  \\
\end{array}\right)^{-1}
```
where the inverse exists because ``T`` is non-degenerate [1], cf. functions
`_cart_to_bary` and `_compute_cart_to_bary_matrix`. Additionally, we have
``\partial_{x_i} λ_j(\boldsymbol{x}) = M_{i+1,j}``, so
```math
\nabla λ_j = M_{2:N, j}.
```
The matrix ``M`` is all we need that depends on ``T`` in order to compute
Bernstein polynomials and their derivatives, it is stored in the field
`cart_to_bary_matrix` of [`BernsteinBasisOnSimplex`](@ref), when ``T`` is not
the reference simplex.

On the reference simplex defined by the vertices `get_vertex_coordinates(SEGMENT / TRI / TET...)`:
```math
\begin{aligned}
v_1     & = (0\ 0\ 0\ ...\ 0), \\
v_2     & = (0\ 1\ 0\ ...\ 0), \\
\vdots  &           \\
v_N     & = (0\ 0\ ...\ 0\ 1),
\end{aligned}
```
the matrix ``M`` is not stored because
```math
λ(\boldsymbol{x}) = \big(1-\sum_{1\leq i\leq D} x_i, x_1, x_2, ..., x_D\big)
\quad\text{and}\quad
\partial_{x_i} λ_j = δ_{i+1,j} - δ_{1j}.
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
B^{D,K}_α(x) = \binom{K}{α} λ^α\qquad\text{for all }α \text{ s.t. } |α|=K,\ \#α=D+1
```
where
- ``|α|=\sum_{1≤ i≤ N} α_i``
- ``\binom{K}{α} = \frac{K!}{α_1 !α_2 !\dotsα_N!}``
- ``λ`` are the barycentric coordinates relative to ``T`` (defined above)

The superscript ``D`` and ``K`` can be omitted because they are always
determined by ``α`` using ``D=\#(α)-1`` and ``K=|α|``.
The set ``\{B_α\}_{α}`` is a basis of
``ℙ^D_K``, implemented by [`BernsteinBasisOnSimplex`](@ref).


### Pascal triangle, pyramid and simplex



### Standard de Casteljau algorithm


### Downward de Casteljau algorithm


### Evaluation of Bernstein polynomials and their derivatives

```math
B_α = \sum_{1≤ i≤ N} λ_i B_{α-e_i} = \sum_{(i,β)\in sub(α)} λ_i B_β
```

```math
\partial_q B_α = K\sum_{|β|=|α|-1} B_β (
δ_{β}^{α-e_q}
-δ_{β}^{α-e_N}
) = K\sum_{(i,β)\in sub(α)}  ( \dots )
```

```math
\partial_t\partial_q B_α = K(K-1)\sum_{|β|=|α|-2} B_β (
δ_{β}^{α-e_t-e_q}
-δ_{β}^{α-e_t-e_N}
-δ_{β}^{α-e_N-e_q}
+δ_{β}^{α-e_N-e_N}
) = K(K-1)\sum_{(i,β)\in sub(α,2)}  ( \dots )
```

```math
```

# Bernstein basis generalization for ``P_r^{(-)}Λ^k`` spaces

The `TODO` basis implements the polynomial bases for the spaces ``ℙ_r^-Λ^k``
and ``ℙ_rΛ^k`` (we write ``P_r^{(-)}Λ^k`` for either one of them) derived in
[2] on simplices of any dimension, for any form degree ``k`` and polynomial
degree ``r``. In this section, we give a summary of the basis definition,
translation of the formulas in the paper from the differential geometry
language to usual differential calculus and the implemented algorithm.

Again, a ``D``-dimensional simplex ``T`` is defined by ``N=D+1`` vertices
``\{v_1, v_2, ..., v_N\}=\{v_i\}_{i\in 1:N}``. We uniquely identify a
``d``-dimensional face ``f`` of ``T`` by the set of the ``d+1`` increasing
indices of its vertices:
```math
f\sim I_f = \{i_0, i_1, ..., i_d\} \qquad\text{such that } 1≤ i_0 < i_1 < ... <i_d≤ N .
```
In particular, ``T\sim \{1:N\}``. We write ``f\subset T`` for any face of
``T``, including ``T`` itself or its vertices.
``T`` has ``\binom{N}{d+1}`` ``d``-dimensional faces, indexed ``\forall\,1≤ i_0
< i_1 < ... < i_d ≤ N``. The dimension of a face ``f`` is ``\#I_f``, and we
write ``{"∀\,\#J=d+1"}`` for all the index sets of the ``d``-dimensional faces of
``T``.

Using Einstein's convention of summation on repeated indices, a degree-``k``
dimension-``D`` form ``ω`` can be written in the canonical Cartesian basis as
``ω = ω_I\,\text{d}x^I``, where the basis is
```math
\big\{ \text{d}x^I = \underset{i\in I}{\bigwedge}\text{d}x^{i}
=\text{d}x^{i_1}\wedge ...\wedge \text{d}x^{i_k} \quad\big|\quad I=\{i_1, ...,
i_k\} \text{ for }1≤ i_1 < ... < i_k ≤ D\big\},
```
``\{ω_I\}_I\in\mathbb{R}^\binom{D}{k}`` is the vector of coefficients of ``ω``,
and ``\{\text{d}x^i\}_{1\leq i\leq D}`` is the canonical covector basis (basis
of ``\text{T}_x T``) such that ``\text{d}x^i(\partial_{x_j})=δ_{ij}``.


#### Geometric decomposition

The main feature of the ``P_r^{(-)}Λ^k`` bases is that each basis polynomial
``ω^{f,α}`` is associated with a face ``f`` of ``T`` (and a domain point
``\boldsymbol{x}_α`` inside ``f``), in the sense that the trace of ``ω^{f,α}``
on another face ``g\subset T`` is zero when ``g`` does not contain ``f``:
```math
f\not\subset g\ \rightarrow\ \text{tr}_g ω^{f,α} = 0, \quad\forall f,g \subseteq T,\ \forall \llbracket α\rrbracket\subseteq I_f, α>0,
```
including any face ``g\neq f`` of dimension less or equal that of ``f``.

These basis polynomials ``ω^{f,α}`` are called bubble functions associated to
``f``, the space they span is called ``\mathring{ℙ}_r^{(-)}Λ^k(T,f)``. There
are no bubble functions of degree ``k`` on faces of dimension ``<k``, so the
spaces ``ℙ_r^{(-)}Λ^k(T)`` admit the geometric decomposition:
```math
ℙ_r^{(-)}Λ^k(T) = \underset{f\subset T}{\oplus}\ \mathring{ℙ}_r^{(-)}Λ^k(f)
= \underset{k≤d≤D}{\oplus}\underset{\quad I_f=0≤ i_0 < ... < i_d ≤ D}{\oplus}\ \mathring{ℙ}_r^{(-)}Λ^k(T,f).
```

#### Bubble functions ``\mathring{ℙ}_r^-Λ^k``

The ``ℙ^-`` type bubble basis polynomials associated to a face ``f\subset T``
are defined by
```math
\mathring{ℙ}_r^-Λ^k(T,f) = \text{span}\big\{B_α φ^J \quad\big|\quad |α|=r-1,\ \#J=k+1,\ ⟦α⟧∪J=I_f,\ α_i=0 \text{ if } i< \text{min}(J) \big\}\newline
```
where ``B_α`` are the scalar Bernstein polynomials implemented by
[`BernsteinBasisOnSimplex`](@ref), and ``φ^J`` are the Whitney forms:
```math
φ^J = \sum_{0≤l≤k} (-1)^{l} λ_l \, \text{d}λ^{J\backslash l} \quad\text{where}\quad
\text{d}λ^{J\backslash l} = \underset{j\in J\backslash \{J_l\} }{\bigwedge}\text{d}λ^{j},
```
``φ^J `` is a ``k``-form of polynomial order ``1``. We now need to express
``\text{d}λ^{J\backslash l}`` in the Cartesian basis ``{\text{d}x^I}``.

In a polytopal simplex ``T`` (flat faces), the 1-forms
``\text{d}λ^j:=\text{d}(λ_j)`` is homogeneous, its coefficients in the
canonical basis are derived by
```math
\text{d}λ^j = (\nabla(λ_j))^♭ = δ_{ki}\partial_{k}λ_j\,\text{d}x^i = \partial_{i}λ_j\,\text{d}x^i = M_{i+1,j}\,\text{d}x^i
```
where ``{}^♭`` is the flat map, the metric ``g_{ki}=δ_{ki}`` is trivial and
``M_{i+1,j}`` are components of the barycentric change of coordinate matrix
``M`` introduced in the Barycentric coordinates section above.

So the exterior products ``\text{d}λ^{J\backslash l}`` are expressed using
determinants of ``k``-minors of ``M`` as follows:
```math
    \text{d}λ^{J\backslash l} = m_I^{J\backslash l}\text{d}x^I
\quad\text{where}\quad m_I^{J\backslash l}
    = \text{det}\big(M_I^{J\backslash l}\big)
    = \text{det}(M_{i+1,j})_{\{i\in I,\, j\in{J\backslash l}\}},
```
and finally the coordinates of ``ω^{α,J}:=B_α φ^J`` in the basis ``\mathrm{d}x^I`` are
```math
ω_I^{α,J} = B_α \sum_{0≤l≤k} (-1)^{l} λ_l \, m_I^{J\backslash l}.
```
There are ``\binom{D}{k}`` coordinates. The ``\binom{D}{k}\binom{N}{k}``
coefficients ``\{m_I^{J}\}_{I,J}`` are constant and are pre-computed from ``M``
at the creation of the basis `TODO`.

Finally, the pseudocode to evaluate our basis of ``ℙ_r^-Λ^k(T)`` at
``\boldsymbol{x}`` is
```
compute λ(x)
compute B_α(x) for all |α|=r-1

for d in k:D
    for #f = k, f≤N
        for α,J in bubble_indices(k,f)
            J' = [J\J1, J\J2, ..., J\Jk]
            for #I = k, I≤D
                s = 0
                for l in 1:k
                    Jl = J'[l]
                    s += (-1)^l * λ_l * m[I][Jl]
                end
                ω^{α,J}_I = B_α * s
            end
            set_value ω^{α,J}
        end
    end
end
```
The loops are unrolled and the indices are pre-computed at compile time using a
`@generated` function barrier in `TODO`.

#### Bubble functions ``\mathring{ℙ}_rΛ^k``

The ``ℙ^-`` type bubble basis polynomials associated to a face ``f\subset T``
are defined by
```math
\mathring{ℙ}_rΛ^k(T,f) = \text{span}\big\{B_α ψ_J^{α,f} \quad\big|\quad |α|=r,\ \#J=k,\ ⟦α⟧∪J=I_f,\ α_i=0 \text{ if } i< \text{min}(I_f \backslash J) \big\},
```
where ``ψ_J^{α,f}`` are defined by
```math
    ψ^{α,J,f} = \underset{j\in J}{\bigwedge} ψ^{α,j,f}
\quad\text{where}\quad
    ψ^{α,j,f} = \mathrm{d}λ^j - \frac{α_j}{|α|}\sum_{i\in I_f}\mathrm{d}λ^i.
```

```math
```

## References

[1] [M.J. Lai & L.L. Schumaker, Spline Functions on Triangulations, Chapter 2 - Bernstein–Bézier Methods for Bivariate Polynomials, pp. 18 - 61.](https://doi.org/10.1017/CBO9780511721588.003)

[2] [D.N. Arnold, R.S. Falk & R. Winther, Geometric decompositions and local bases for spaces of finite element differential forms, Computer Methods in Applied Mechanics and Engineering](https://doi.org/10.1016/j.cma.2008.12.017)


