
# Bernstein bases algorithms

### Barycentric coordinates

A ``D``-dimensional simplex ``T`` is defined by ``N=D+1`` vertices ``\{v_1,
v_2, ..., v_N\}=\{v_i\}_{i\in 1:N}``. The barycentric coordinates
``λ(\boldsymbol{x})=\{λ_j(\boldsymbol{x})\}_{1≤ j≤ N}`` are uniquely
defined by:
```math
\boldsymbol{x} = \sum_{1≤ j≤ N} λ_j(\boldsymbol{x})v_j \quad\text{and}\quad
\sum_{1≤ j≤ N} λ_j(\boldsymbol{x}) = 1,
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
``\partial_{x_i} λ_j(\boldsymbol{x}) = M_{j,i+1}``, so
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
λ(\boldsymbol{x}) = \big(1-\sum_{1≤ i≤ D} x_i, x_1, x_2, ..., x_D\big)
\quad\text{and}\quad
\partial_{x_i} λ_j = δ_{i+1,j} - δ_{1j} = M_{j,i+1}.
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

# Bernstein basis generalization for ``ℙ_r^{(-)}Λ^k`` spaces

The `TODO` basis implements the polynomial bases for the spaces ``ℙ_r^-Λ^k(T^D)``
and ``ℙ_rΛ^k(T^D)`` (we write ``ℙ_r^{(-)}Λ^k`` for either one of them) derived in
[2] on simplices of any dimension, for any form degree ``k`` and polynomial
degree ``r``. In this section, we give a summary of the basis definition,
translation of the formulas in the paper from the differential geometry
language to usual differential calculus and the implemented algorithm.

#### Face and form coefficients indexing

Again, a ``D``-dimensional simplex ``T`` is defined by ``N=D+1`` vertices
``\{v_1, v_2, ..., v_N\}=\{v_i\}_{i\in 1:N}``. We uniquely identify a
``d``-dimensional face ``f`` of ``T`` by the set ``F`` of the ``d+1``
increasing indices of its vertices:
```math
f\sim F = \{F_0, F_1, ..., F_d\} \qquad\text{such that } 1≤ F_0 < F_1 < ... <F_d≤ N .
```
In particular, ``T\sim \{1:N\}``. We write ``f\subset T`` for any face of
``T``, including ``T`` itself or its vertices.
``T`` has ``\binom{N}{d+1}`` ``d``-dimensional faces, indexed ``\forall\,1≤ F_0
< F_1 < ... < F_d ≤ N``. The dimension of a face ``f`` is ``\#F``, and we
write ``{"∀\,\#J=d+1"}`` for all the increasing index sets of the
``d``-dimensional faces of ``T``.

Using Einstein's convention of summation on repeated indices, a degree-``k``
dimension-``D`` form ``ω`` can be written in the canonical Cartesian basis as
``ω = ω_I\,\text{d}x^I``, where the basis is
```math
\big\{ \text{d}x^I = \underset{i\in I}{\bigwedge}\text{d}x^{i}
=\text{d}x^{I_1}\wedge ...\wedge \text{d}x^{I_k} \quad\big|\quad I=\{I_1, ...,
I_k\} \text{ for }1≤ I_1 < ... < I_k ≤ D\big\},
```
``\{ω_I\}_I\in\mathbb{R}^\binom{D}{k}`` is the vector of coefficients of ``ω``,
and ``\{\text{d}x^i\}_{1≤ i≤ D}`` is the canonical covector basis (basis
of ``\text{T}_x T``) such that ``\text{d}x^i(\partial_{x_j})=δ_{ij}``.

TODO choose linear ordering

#### Geometric decomposition

The main feature of the ``ℙ_r^{(-)}Λ^k`` bases is that each basis polynomial
``ω^{f,α}`` is associated with a face ``f`` of ``T`` (and a domain point
``\boldsymbol{x}_α`` inside ``f``), in the sense that the trace of ``ω^{f,α}``
on another face ``g\subset T`` is zero when ``g`` does not contain ``f``:
```math
f\not\subset g\ \rightarrow\ \text{tr}_g ω^{f,α} = 0, \quad\forall f,g \subseteq T,\ \forall \llbracket α\rrbracket\subseteq F, α>0,
```
including any face ``g\neq f`` of dimension less or equal that of ``f``.

These basis polynomials ``ω^{f,α}`` are called bubble functions associated to
``f``, the space they span is called ``\mathring{ℙ}_r^{(-)}Λ^k(T,f)``. There
are no bubble functions of degree ``k`` on faces of dimension ``<k``, so the
spaces ``ℙ_r^{(-)}Λ^k(T)`` admit the geometric decomposition:
```math
ℙ_r^{(-)}Λ^k(T) = \underset{f\subset T}{\oplus}\ \mathring{ℙ}_r^{(-)}Λ^k(f)
= \underset{k≤d≤D}{\oplus}\underset{\quad F=0≤ i_0 < ... < i_d ≤ D}{\oplus}\ \mathring{ℙ}_r^{(-)}Λ^k(T,f).
```

#### Bubble functions ``\mathring{ℙ}_r^-Λ^k``

The ``ℙ^-`` type bubble basis polynomials associated to a face ``f\subset T``
are defined by
```math
\mathring{ℙ}_r^-Λ^k(T,f) = \text{span}\big\{ \bar{ω}_α^J = B_α φ^J \quad\big|\quad |α|=r-1,\ \#J=k+1,\ ⟦α⟧∪J=F,\ α_i=0 \text{ if } i< \text{min}(J) \big\}\newline
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
\text{d}λ^j = (\nabla(λ_j))^♭ = δ_{ki}\partial_{k}λ_j\,\text{d}x^i =
\partial_{i}λ_j\,\text{d}x^i = M_{j,i+1}\,\text{d}x^i
```
where ``{}^♭`` is the flat map, the metric ``g_{ki}=δ_{ki}`` is trivial and
``M_{j,i+1}`` are components of the barycentric change of coordinate matrix
``M`` introduced in the Barycentric coordinates section above.

So the exterior products ``\text{d}λ^{J\backslash l}`` are expressed using
determinants of ``k``-minors of ``M`` as follows:
```math
    \text{d}λ^{J\backslash l} = m_I^{J\backslash l}\text{d}x^I
\quad\text{where}\quad m_I^J
    = \text{det}\big( (\partial_{I_i}λ_{J_j})_{1\leq i,j\leq k} \big)
    = \text{det}\big( (M_{J_j,I_i+1})_{1\leq i,j\leq k} \big),
```
and we obtain the coordinates of ``\bar{ω}_α^{J}=B_α φ^J`` in the basis
``\mathrm{d}x^I`` are
```math
\bar{ω}_{α,I}^{J} = B_α \sum_{0≤l≤k} (-1)^{l} λ_l \, m_I^{J\backslash l}.
```
There are ``\binom{D}{k}`` coordinates. The ``\binom{D}{k}\binom{N}{k}``
coefficients ``\{m_I^{J}\}_{I,J}`` are constant in ``T`` and are pre-computed
from ``M`` at the creation of the basis `TODO`.

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
                ω_αJ[I] = B_α * s
            end
            set_value(result, ω_αJ})
        end
    end
end
```
The loops are unrolled and the indices are pre-computed at compile time using a
`@generated` function barrier in `TODO`.

#### Bubble functions ``\mathring{ℙ}_rΛ^k``

The ``ℙ`` type bubble basis polynomials associated to a face ``f\subset T``
are defined by
```math
\mathring{ℙ}_rΛ^k(T,f) = \text{span}\big\{ ω_{α,f}^J=B_α Ψ_{α,f}^J
\quad\big|\quad |α|=r,\ \#J=k,\ ⟦α⟧∪J=F,\ α_i=0 \text{ if } i< \text{min}(F
\backslash J) \big\},
```
where ``Ψ_{α,f}^J`` are defined by
```math
    Ψ_{α,f}^J = \underset{j\in J}{\bigwedge} Ψ_{α,f}^j
\quad\text{where}\quad
 Ψ_{α,f}^j = \mathrm{d}λ^j - \frac{α_j}{|α|}\sum_{l\in F}\mathrm{d}λ^l.
```
Again, we need their coordinates in the Cartesian basis ``\mathrm{d}x^I``:
```math
Ψ_{α,f}^j = M_{j,i+1}\mathrm{d}x^i - \frac{α_j}{|α|}\sum_{l\in F}M_{l,i+1}\mathrm{d}x^i
= \big(M_{j,i+1} - \frac{α_j}{|α|}\sum_{l\in F}M_{l,i+1}\big)\mathrm{d}x^i
```
so
```math
Ψ_{α,f}^j = ψ_{α,f,i}^j \mathrm{d}x^i
\quad\text{where}\quad
ψ_{α,f,i}^j = M_{j,i+1} - \frac{α_j}{|α|}\sum_{l\in F}M_{l,i+1}
```
and

```math
Ψ_{α,f}^J = ψ_{α,f,I}^J \mathrm{d}x^I
\quad\text{where}\quad
ψ_{α,f,I}^J = \text{det}\big( (ψ_{α,f,i}^j)_{i\in I,\,j\in J} \big).
```

Finally, the ``\binom{D}{k}`` coordinates of ``ω_{α,f}^{J}=B_α Ψ_{α,f}^J`` in the
basis ``\mathrm{d}x^I`` are
```math
ω_{α,f,I}^{J} = B_α\, ψ_{α,f,I}^J,
```
where the ``\binom{D+r}{k+r}\binom{r+k}{k}\binom{D}{k}
=\mathrm{dim}(ℙ_rΛ^k(T^D))\times\# (\{\mathrm{d}x^I\}_I)`` coefficients
``ψ_{α,f,I}^J`` depend only on ``T`` and are pre-computed at the construction
of the basis `TODO`.

The pseudocode to evaluate the basis can be written in simpler way:
```
compute λ(x)
compute B_α(x) for all |α|=r

for (Ψ_αfJ, α) in Ψ_α_couples
    ω_αfJ = B_α * Ψ_αfJ
    set_value(result, ω_αfJ)
end
```

#### Exterior derivatives of the basis forms


#### Gradient and Hessian of the coefficient vectors

For retro-compatibility with the rest of Gridap, let us derive the formula for
the gradient and hessian of the form coefficient vectors
``\{\bar{ω}_{α,f,I}^{J}\}_I`` and ``\{ω_{α,f,I}^{J}\}_I``. We will express them
in function of the scalar Bernstein polynomial derivatives already implemented
by `BernsteinBasisOnSimplex`.

##### Coefficient vector ``\{\bar{ω}_{α,f,I}^{J}\}_I``

Recall ``\bar{ω}_{α,I}^{J} = B_α \sum_{0≤l≤k} (-1)^{l} λ_l \, m_I^{J\backslash
l}``, the derivatives are not immediate to compute because both ``B_α`` and
``λ_l`` depend on ``\boldsymbol{x}``, let us first use ``B_α λ_l = \frac{α_l +
1}{|α|+1}B_{α+e_l}`` to write the coefficients in Bernstein form as follows
```math
\bar{ω}_{α,I}^{J} = B_α \sum_{0≤l≤k} (-1)^{l} λ_l \, m_I^{J\backslash l} =
\frac{1}{r}\sum_{0≤l≤k} (-1)^{l} (α_l +1)  m_i^{j\backslash l}\, b_{α+e_l},
```
where ``|α|+1`` was replaced with ``r``, the polynomial degree of
``\bar{ω}_{α,I}^{J}``. As a consequence, for any Cartesian coordinate indices
``1\leq p,q\leq D``, we get
```math
\partial_q \bar{ω}_{α,I}^{J} = \frac{1}{r}\sum_{0≤l≤k} (-1)^{l} (α_l +1)  m_i^{j\backslash l}\, \partial_q B_{α+e_l},\\
\partial_t\partial_q \bar{ω}_{α,I}^{J} = \frac{1}{r}\sum_{0≤l≤k} (-1)^{l} (α_l +1)  m_i^{j\backslash l}\, \partial_t\partial_q B_{α+e_l}.
```

##### Coefficient vector ``\{ω_{α,f,I}^{J}\}_I``

Recall ``ω_{α,f,I}^{J} = B_α\, ψ_{α,f,I}^J``, the derivatives are easy to
compute because only ``B_α`` depends on ``\boldsymbol{x}``, leading to
```math
\partial_q\, ω_{α,f,I}^{J} = ψ_{α,f,I}^J\; \partial_q B_α,\\
\partial_t\partial_q\, ω_{α,f,I}^{J} = ψ_{α,f,I}^J\; \partial_t\partial_q B_α.
```

#### Optimizations for the reference simplex

In the reference simplex ``\hat{T}``, the vertices and thus coefficients of
``M`` are known at compile time, so the coefficients ``m_I^J`` and
``ψ_{α,f,I}^J`` in ``\hat{T}``, denoted by ``\hat{m}_I^J`` and
``\hat{ψ}_{α,f,I}^J`` respectively, can be hard-coded at compile time in the
`@generated` functions to avoid storing them in the basis and accessing them at
runtime. Let us derive the formulas for them.

##### Coefficients ``\hat{m}_I^J``

It was shown in the Barycentric coordinates section above that
``M_{j,i+1} = δ_{i+1,j} - δ_{1j}``. Let ``\#I=k`` and ``\#J=k``. We need
to compute the determinant of the matrix
```math
\hat{M}_{IJ}=(δ_{I_i+1,J_j}-δ_{1,J_j})_{1\leq i,j\leq k}.
```
Let us define:
- ``s=δ_1^{J_1}``, that indicates if ``\hat{M}_{IJ}`` contains a column of ``-1``,
- ``p = \text{min } \{j\,|\, I_j+s \notin J\}`` where ``\text{min}\,\emptyset=0``, the index of the first row of ``\hat{M}_{IJ}`` containing no ``1``. ``p=0`` if and only if ``\hat{M}_{IJ}`` is the identity matrix.
- ``n`` be the number of columns of zeros of ``\hat{M}_{IJ}``, defined by ``n =\# \{\ i\ |\ J_i-s\notin I\}``.
Then it can be shown that ``k-n`` is the rank of ``\hat{M}_{IJ}``, and that
```math
\hat{m}_I^J = \mathrm{det}(\hat{M}_{IJ}) = (-1)^pδ_{n,0},
```
so in ``\hat{T}``, there is
```math
\bar{ω}_{α,I}^{J} = B_α \sum_{0≤l≤k} (-1)^{l} λ_l \, \hat{m}_I^{J\backslash l}.
```

##### Coefficients ``\hat{ψ}_{α,f,I}^J``
The expression of ``ψ_{α,f,i}^j`` in ``\hat{T}`` is
```math
\hat{ψ}_{α,f,i}^j = M_{j,i+1} - \frac{α_j}{|α|}\sum_{l\in F}δ_{i+1,l} - δ_{1l}
= M_{j,i+1} + \big( δ_{1,F_1}-\sum_{l\in F}δ_{i+1,l} \big)\frac{α_j}{|α|}
```
leading to
```math
\hat{ψ}_{α,f,I}^J  = \mathrm{det}\Big(\hat{M}_{IJ} + u\,v^{\mathrm{T}}\Big)
\quad\text{where}\quad
u^i = δ_{1,F_1}-\sum_{l\in F}δ_{I_i+1,l}, \qquad v^j = \frac{α_{J_j}}{|α|}.
```
We can use the following matrix determinant lemma:
```math
\mathrm{det}(\hat{M}_{IJ} + uv^T) =
\mathrm{det}(\hat{M}_{IJ}) + v^T\mathrm{adj}(\hat{M}_{IJ})u.
```
The determinant ``\mathrm{det}(\hat{M}_{IJ})=\hat{m}_I^{J}`` was computed
above, but ``\mathrm{adj}(\hat{M}_{IJ})``, the transpose of the cofactor matrix
of ``\hat{M}_{IJ}``, is also needed. Let ``s=δ_1^{J_1}``, ``n`` and ``p`` be
defined as above, and additionally define
- ``q = \text{min } \{j\,|\,j>p,\ I_j+s \notin J\}``, the index of the second row of ``\hat{M}_{IJ}`` containing no ``1`` (``q=0`` if there isn't any),
- ``m = \text{min } \{i\,|\,i>s,\ J_i-s \notin I\}``, the index of the first column of ``\hat{M}_{IJ}`` containing only zeros (``m=0`` if there isn't any).

Then the following table gives the required information to apply the matrix
determinant lemma and formulas for ``\hat{ψ}_{α,f,I}^J``
```math
\begin{array}{|c|c|c|c|c|}
\hline
s  & n & \mathrm{rank}\hat{M}_{IJ} & \mathrm{adj}\hat{M}_{IJ} & \hat{ψ}_{α,f,I}^J \\
\hline
\hline
0   & 0 & k   & δ_{ij}                         & 1 + u \cdot v\\
\hline
0   & 1 & k-1 & (-1)^{m+p}δ_i^m δ^p_j          & (-1)^{m+p}u^m v^p \\
\hline
1   & 0 & k   & (-1)^p(δ_{J_i,I_j+1}-δ_{p,J_j})& (-1)^p(1-u^p|v|+\underset{1\leq l<p}{\sum}v^{l+1}u^l + \underset{p<l\leq k}{\sum}v^{l}u^l) \\
\hline\hspace{1mm}
1   & 1 & k-1 & (-1)^{m+p+q}δ_i^m(δ^q_j-δ^p_j) & (-1)^{m+p+q}u^m(v^q-v^p) \\
\hline
0/1 & \geq 2 & \leq k-2 & 0 & 0 \\
\hline
\end{array}
```
In this table, ``m``, ``p`` and ``q`` depend on ``I`` and ``J``, ``u``
depends on ``F`` and ``I``, and ``v`` depends on ``α`` and ``J``.

```math
```


#### Trace of the basis forms on the simplex faces

TODO

#### Hodge operator of the basis forms

TODO

##### Useful lemmas TODO

##### Lemma 1
For ``I`` and ``J`` two increasing sets of indices of same size ``k``,
```math
\mathrm{det}\big( (δ_{I_i, J_j})_{1\leq i,j\leq k} \big) = δ_I^J \ \text{(}= \underset{1\leq i\leq k}{\Pi} δ_{I_i}^{J_i} \text{)}
```

Proof:
Suppose ``I_1\neq J_1``, because both sets are strictly increasing, either
``I_1 < J_1 < J_j`` and the first row of the matrix contains only zeros, or
``I_i > I_1 > J_1`` and the first column contains only zeros. So in both cases,
the determinant would be zero, otherwise one can do a Laplace expansion along
the first row, giving
```math
\mathrm{det}\big( (δ_{I_i, J_j})_{1\leq i,j\leq k} \big) = δ_{I_1}^{J_1} \mathrm{det}\big((δ_{I_i, J_j})_{2\leq i,j\leq k}\big).
```
The same argument is repeated recursively to obtain the result.

##### Lemma X

Let us determine the rank of ``A`` from ``I`` and ``J``, leading to 6
possibilities for the expression of the adjoint depending on weather the rank
of ``A`` is ``k``, ``k-1`` or ``<k``, and weather ``J_1`` is ``1``.
If ``\mathrm{rank}(A)\leq k-2``, both ``\mathrm{det}(A)`` and
``\mathrm{adj}(A)`` vanish.

## References

[1] [M.J. Lai & L.L. Schumaker, Spline Functions on Triangulations, Chapter 2 - Bernstein–Bézier Methods for Bivariate Polynomials, pp. 18 - 61.](https://doi.org/10.1017/CBO9780511721588.003)

[2] [D.N. Arnold, R.S. Falk & R. Winther, Geometric decompositions and local bases for spaces of finite element differential forms, Computer Methods in Applied Mechanics and Engineering](https://doi.org/10.1016/j.cma.2008.12.017)


