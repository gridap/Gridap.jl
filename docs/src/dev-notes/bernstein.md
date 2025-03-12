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
v_2     & = (0\ 1\ ⋯\ 0), \\
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

# Bernstein basis generalization for ``ℙ_r^{(-)}Λ^k`` spaces

The `TODO` basis implements the polynomial bases for the spaces ``ℙ_r^-Λ^k(T^D)``
and ``ℙ_rΛ^k(T^D)`` (we write ``ℙ_r^{(-)}Λ^k`` for either one of them) derived in
[2] on simplices of any dimension, for any form degree ``k`` and polynomial
degree ``r``. In this section, we give a summary of the basis definition,
translation of the formulas in the paper from the differential geometry
language to usual differential calculus and the implemented algorithm.

Again, a ``D``-dimensional simplex ``T`` is defined by ``N=D+1`` vertices
``\{v_1, v_2, ..., v_N\}=\{v_i\}_{i\in 1:N}``. We uniquely identify a
``d``-dimensional face ``f`` of ``T`` by the set ``F`` of the ``d+1``
increasing indices of its vertices:
```math
f\sim F = \{f_0, f_1, ..., f_d\} \qquad\text{such that } 1≤ f_0 < f_1 < ... <f_d≤ N .
```
In particular, ``T\sim \{1:N\}``. We write ``f\subset T`` for any face of
``T``, including ``T`` itself or its vertices.
``T`` has ``\binom{N}{d+1}`` ``d``-dimensional faces, indexed ``\forall\,1≤ f_0
< f_1 < ... < f_d ≤ N``. The dimension of a face ``f`` is ``\#F``, and we
write ``{"∀\,\#J=d+1"}`` for all the increasing index sets of the
``d``-dimensional faces of ``T``.

Using Einstein's convention of summation on repeated indices, a degree-``k``
dimension-``D`` form ``ω`` can be written in the canonical Cartesian basis as
``ω = ω_I\,\text{d}x^I``, where the basis is
```math
\big\{ \text{d}x^I = \underset{i\in I}{\bigwedge}\text{d}x^{i}
=\text{d}x^{i_1}\wedge ...\wedge \text{d}x^{i_k} \quad\big|\quad I=\{i_1, ...,
i_k\} \text{ for }1≤ i_1 < ... < i_k ≤ D\big\},
```
``\{ω_I\}_I\in\mathbb{R}^\binom{D}{k}`` is the vector of coefficients of ``ω``,
and ``\{\text{d}x^i\}_{1≤ i≤ D}`` is the canonical covector basis (basis
of ``\text{T}_x T``) such that ``\text{d}x^i(\partial_{x_j})=δ_{ij}``.


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
ω_{α,f,I}^{J} = B_α ψ_{α,f,I}^J,
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

#### Optimizations for the reference simplex

In the reference simplex ``\hat{T}``, the vertices and thus coefficients of
``M`` are known at compile time, so the coefficients ``m_I^J`` and
``ψ_{α,f,I}^J`` in ``\hat{T}``, denoted by ``\hat{m}_I^J`` and
``\hat{ψ}_{α,f,I}^J`` respectively, can be hard-coded at compile time in the
`@generated` functions to avoid storing them in the basis and accessing them at
runtime. Let us derive the formulas for them.

##### Coefficients ``\hat{m}_I^J``

It was shown in the Barycentric coordinates section above that ``M_{j,i+1} =
δ_{i+1,j} - δ_{1j}``. Let ``\#I=k`` and ``\#J=k``, and let
``m = \underset{1\leq i\leq k-1}{\text{argfirst }}(I_i+1\neq J_i)`` or
``m=k`` if no different index is found. Then it can be shown that
```math
\hat{m}_I^J = \mathrm{det}\big((δ_{I_i+1,J_j}-δ_{1,J_j})_{1\leq i,j\leq k}\big) =

\left\{\begin{array}{cl}
    (-1)^m δ_{I\backslash m \;+ 1}^{J\backslash 1}& \text{ if } J_1 = 1\\
δ_{I+1}^{J} &\text{ else }
\end{array}\right.
,
```
where ``I+1=\{I_1+1, ...,I_k+1\}`` and ``δ_I^J = \underset{1\leq l\leq
k}{\Pi}δ_{I_l}^{J_l}``. So in ``\hat{T}``, there is
```math
\bar{ω}_{α,I}^{J} = B_α \sum_{0≤l≤k} (-1)^{l} λ_l \, \hat{m}_I^{J\backslash l}.
```

##### Coefficients ``\hat{ψ}_{α,f,I}^J``
```math
\hat{ψ}_{α,f,i}^j = δ_{i+1,j} - δ_{1j} - \frac{α_j}{|α|}\sum_{l\in F}δ_{i+1,l} - δ_{1l}
= δ_{i+1,j} - δ_{1j} + \big( δ_{1,F_1}-\sum_{l\in F}δ_{i+1,l} \big)\frac{α_j}{|α|}
```

```math
\hat{ψ}_{α,f,I}^J = \mathrm{det}\Big( \Big(
    δ_{I_i+1,J_j}-δ_{1,J_j} - \big( δ_{1,F_1}-\sum_{l\in F}δ_{I_i+1,l} \big)\frac{α_{J_j}}{|α|}
\Big)_{1\leq i,j\leq k} \Big)
```

We can use the matrix determinant lemma ``\mathrm{det}(A + uv^T) = \mathrm{det}(A) + v^T\mathrm{adj}(A)u`` with
```math
A_{ij}=(δ_{I_i+1,J_j}-δ_{1,J_j})_{1\leq i,j\leq k},\qquad u^i = δ_{1,F_1}-\sum_{l\in F}δ_{I_i+1,l}, \qquad v^j = \frac{α_{J_j}}{|α|}.
```
The determinant ``\mathrm{det}(A)=\hat{m}_I^{J}`` was computed above, but
``\mathrm{adj}(A)``, the transpose of the cofactor matrix of ``A``, is also
needed. Let ``s=δ_1^{J_1}`` and ``m,p,q`` defined by
```math
m = \text{min } \{i\,|\, J_i-s \notin I,\,i>s\},\ \quad
p = \text{min } \{j\,|\, I_j+s \notin J\}  \quad\text{and}\quad
q = \text{min } \{j\,|\, I_j+s \notin J,\,j>p\},
```
where ``\text{min}\,\emptyset=0``, and let ``n`` be the number
of columns of zeros of ``A``, defined by ``n =\# \{\ i\ |\ J_i-s\notin I\}``.
Then it can be shown that ``\text{rank}(A)=k-n+s``, and the following table
gives the required information to apply the matrix determinant lemma and
formulas for ``\hat{ψ}_{α,f,I}^J``
```math
\begin{array}{|c|c|c|c|c|}
\hline
s = δ_1^{J_1} & n & \mathrm{rank}(A) & \mathrm{adj}(A) & \hat{ψ}_{α,f,I}^J \\
\hline
0   & 0       & k       & δ_{ij}                         & 1 + u \cdot v\\
\hline
0   & 1       & k-1     & (-1)^{m+p}δ_i^m δ^p_j          & (-1)^{m+p}u^m v^p \\
\hline
1   & 1       & k       &                          & (-1)^p(1-u^p|v|+\underset{1\leq l<p}{\sum}v^{l+1}u^l + \underset{p<l\leq k}{\sum}v^{l}u^l) \\
\hline\hspace{1mm}
1   & 2       & k-1     & (-1)^{m+p+q}δ_i^m(δ^q_j-δ^p_j) & (-1)^{m+p+q}u^m(v^q-v^p) \\
\hline
0/1 & 2+s\leq & k-2\geq & 0 & 0 \\
\hline
\end{array}
```
In this table, ``m``, ``p`` and ``q`` depend on ``I`` and ``J``, ``u``
depends on ``F`` and ``I``, and ``v`` depends on ``α`` and ``J``.

```math
```

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


