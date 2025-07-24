```@meta
CurrentModule = Gridap.Polynomials
```

# Bernstein bases algorithms

## Barycentric coordinates

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
[`_cart_to_bary`](@ref) and [`_compute_cart_to_bary_matrix`](@ref). Additionally,
we have ``∂_{x_i} λ_j(\bm{x}) = M_{j,i+1}``, so
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

## Bernstein polynomials definition

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

## Bernstein basis generalization for ``ℙΛ`` spaces

The [`PmLambdaBasis`](@ref) and [`PLambdaBasis`](@ref) bases respectively
implement the polynomial bases for the spaces ``ℙ_r^-Λ^k(T^D)`` and
``ℙ_rΛ^k(T^D)`` (we write ``ℙ_r^{(-)}Λ^k`` for either one of them) derived in
[2] on simplices of any dimension, for any form degree ``k`` and polynomial
degree ``r``. These spaces include and generalize several standard FE
polynomial spaces, see the Periodic Table of the Finite Elements [3].

The following notes explain the implementation in detail. For the moment, only
the space with form order ``{k = 0,1,D-1}`` and ``D`` are available, because the
forms are translated into their [vector calculus proxy](@ref "Translation between forms and vectors").

#### Face and form coefficients indexing

Again, a ``D``-dimensional simplex ``T`` is defined by ``N=D+1`` vertices
``\{v_1, v_2, ..., v_N\}=\{v_i\}_{i∈1:N}``. We uniquely identify a
``d``-dimensional face ``F`` of ``T`` by the set of the ``d+1``
increasing indices of its vertices:
```math
F = \{F_1, F_2, ..., F_{d+1}\}
\qquad\text{such that } 1≤ F_1 < F_2 < ... <F_{d+1}≤ N .
```
In particular, ``T\sim \{1:N\}``. We write ``F⊆ T`` for any face of ``T``,
including ``T`` itself or its vertices. ``T`` has ``\binom{N}{d+1}``
``d``-dimensional faces, indexed ``\forall\,1≤ F_1 < F_2 < ... < F_{d+1} ≤ N``.
The dimension of a face ``F`` is ``\#F\;`` (`length(F)`), and we write
``{"∀\,\#J=d+1"}`` for all the increasing index sets of the ``d``-dimensional
faces of ``T``. We will sometimes write ``_{J(i)}`` instead of ``_{J_i}`` for
readability purpose when it appears as a subscript.

Using Einstein's convention of summation on repeated indices, a degree-``k``
dimension-``D`` form ``ω`` can be written in the canonical Cartesian basis as
``ω = ω_I\,\text{d}x^I``, where the basis is
```math
\big\{ \text{d}x^I = \underset{i∈I}{⋀}\text{d}x^{i}
=\text{d}x^{I_1}∧ ...∧ \text{d}x^{I_k} \quad\big|\quad I=\{I_1, ...,
I_k\} \text{ for }1≤ I_1 < ... < I_k ≤ D\big\},
```
``\{ω_I\}_I∈\mathbb{R}^\binom{D}{k}`` is the vector of coefficients of ``ω``,
and ``\{\text{d}x^i\}_{1≤ i≤ D}`` is the canonical covector basis (basis
of ``\text{T}_x T``) such that ``\text{d}x^i(∂_{x_j})=δ_{ij}``.

These sets of indices ``I,J,F`` are ``k``-combinations of ``{1:D/N}``, stored
in `Vector{Int}`. A generator [`_sorted_combinations`](@ref) returns a vector
containing all the ``D``-dimensional ``k``-combinations, and
[`_combination_index`](@ref) can be used to compute the index of a combination
in this vector. This `k`-combinations ordering defines the indices of form
components ``ω_I`` in numerical collections. The order is independent of the
dimension.

#### Translation between forms and vectors

By default, the polynomial forms of order ``k = 0,1,D-1`` and ``D`` are
translated into their equivalents in the standard vector calculus framework
(assuming the simplex is Euclidean).

| ``k``    | Form value                               | Vector proxy value          | Proxy value type |
| :------- | :--------------------------------------- | :-------------------------- | :--------------- |
| ``0  ``  | ``ω ∈ \mathbb{R}``                       | ``ω♯ = ω``                  | `T`              |
| ``1  ``  | ``ω = ω_i \mathrm{d}x^i``                | ``ω♯ = \sum_i ω_i \bm{e}_i``|`VectorValue{D,T}`|
| ``D-1≥1``| ``ω=ω_{I}\mathrm{d}x^{I}`` | ``(⋆ω)♯ =\ \underset{I=\{1:D\}\backslash\{i\}}{\sum_i} (-1)^{i+1} ω_{I} \boldsymbol{e}_i`` |`VectorValue{D,T}`|
| ``D  ``  | ``ω=ω_{\{1:D\}}\mathrm{d}x^{\{1:D\}} ``  | ``(⋆ω)♯=ω_{\{1:D\}}``       |`T`               |

This change of coordinate is implemented by [`_basis_forms_components`](@ref),
the indices of a basis `b::P(m)LambdaBasis` are stored in `b._indices.components`.
For ``{D=2}`` and ``{k=1}``, the default proxy is ``ω♯``. The user may choose the
``(⋆ω)♯`` proxy (for div-conforming spaces) using the kwarg `rotate_90=true`.

#### Geometric decomposition

The main feature of the `P(m)LambdaBasis` bases is that each basis polynomial
``ω^{α,J}`` is associated with a face ``F`` of ``T`` via ``{F=⟦α⟧∪J}`` with
``J`` a face of ``F`` and ``α`` a Bernstein index whose associated domain point
``\boldsymbol{x}_α`` is geometrically inside ``F``. Importantly, the trace of
``ω^{α,J}`` on another face ``G⊆ T`` is zero when ``G`` does not contain
``F``:
```math
F\not⊆ G\ \rightarrow\ \text{tr}_G\, ω^{α,J} = 0, \quad\forall F,G ⊆ T,\ \forall α,J \text{ s.t. }\llbracket α\rrbracket\cup J = F,
```
including any face ``G\neq F`` of dimension less or equal that of ``F``.

These basis polynomials ``ω^{α,J}`` are called bubble functions associated to
``F``, the space they span is called ``\mathring{ℙ}_r^{(-)}Λ^k(T,F)``. There
are no bubble functions of degree ``k`` on faces of dimension ``<k``, so the
spaces ``ℙ_r^{(-)}Λ^k(T)`` admit the geometric decomposition:
```math
ℙ_r^{(-)}Λ^k(T) = \underset{F⊆ T}{\oplus}\ \mathring{ℙ}_r^{(-)}Λ^k(F)
= \underset{k≤d≤D}{\oplus}\underset{\quad F=1≤ F_1 < ... < F_{d+1} ≤ N}{\oplus}\ \mathring{ℙ}_r^{(-)}Λ^k(T,F).
```

#### Bubble functions ``\mathring{ℙ}_r^-Λ^k``

The ``ℙ^-`` type bubble basis polynomials associated to a face ``F⊆T``
defined by [2, Th. 6.1-4] are
```math
\mathring{ℙ}_r^-Λ^k(T,F) = \text{span}\big\{ ω̄^{α,J} =
B_α φ^J \ \big| \ α∈\mathcal{I}_{r-1}^D,\ \#J=k\!+\!1,\ ⟦α⟧∪J=F,\ α_i=0 \text{ if } i< \text{min}(J) \big\}
```
where ``B_α`` are the scalar Bernstein polynomials implemented by
[`BernsteinBasisOnSimplex`](@ref), and ``φ^J`` [2, Eq. (6.3)] are the Whitney
forms:
```math
φ^J = \sum_{1≤l≤k+1} (-1)^{l+1} λ_{J(l)} \, \text{d}λ^{J\backslash l} \quad\text{where}\quad
\text{d}λ^{J\backslash l} = \underset{j∈J\backslash \{J_l\} }{⋀}\text{d}λ^{j},
```
``φ^J `` is a ``k``-form of polynomial order ``1``.

Given ``k,r`` and ``D``, the function [`PmΛ_bubbles(r,k,D)`](@ref PmΛ_bubbles)
computes, for each ``d``-face``F`` "owning" a bubble space, the indices
necessary to compute its bubble polynomials. `PmΛ_bubbles` is used as follows:
```julia
for (F, bubble_functions) in PΛ_bubbles(r,k,D)  # d = length(F)
    for (w, α, α_id, J, sub_J_ids, sup_α_ids) in bubble_functions
        # do stuff for ω̄^{α,J}
    end
end
```
where
- `w` is the index of ``ω̄^{α,J}`` in the whole `PmLambdaBasis`,
- `α` is a `Vector{Int}`,
- `α_id` is [`bernstein_term_id(α)`](@ref bernstein_term_id), the index of `Bα` in the scalar [`BernsteinBasisOnSimplex`](@ref),
- `J` is a `Vector{Int}`,
- `sub_J_ids` is a `::Vector{Int}` are the [`_combination_index`](@ref) of each ``J\backslash \{J(l)\}`` for ``1\leq l\leq \#J``,
- `sup_α_ids` is a `::Vector{Int}` are the [`bernstein_term_id`](@ref) of each ``α+e_i`` for ``1\leq i\leq \#α``.

The implementation is flexible enough to select a subset of the bubble spaces,
the bubbles of a `b::PmLambdaBasis` are obtained via [`get_bubbles(b)`](@ref
get_bubbles) (do NOT modify them).

We now need to express ``\text{d}λ^{J\backslash l}`` in the Cartesian basis
``{\text{d}x^I}``. In a polytopal simplex ``T`` (flat faces), the 1-forms
``\text{d}λ^j:=\text{d}(λ_j)`` is homogeneous, its coefficients in the
canonical basis are derived by
```math
\text{d}λ^j = (\nabla(λ_j))^♭ = δ_{ki}∂_{k}λ_j\,\text{d}x^i =
∂_{i}λ_j\,\text{d}x^i = M_{j,i+1}\,\text{d}x^i
```
where ``{}^♭`` is the flat map, the metric ``g_{ki}=δ_{ki}`` is trivial and
``M_{j,i+1}`` are components of the barycentric change of coordinate matrix
``M`` introduced in the Barycentric coordinates section above.

So the exterior products ``\text{d}λ^{J\backslash l}`` are expressed using
the ``k``-minors ``m_I^{J\backslash l}`` of ``M^\intercal`` as follows:
```math
\text{d}λ^{J\backslash l} = m_I^{J\backslash l}\text{d}x^I
\quad\text{where}\quad m_I^J
= \text{det}\big( (∂_{I(i)}λ_{J(j)})_{1≤ i,j≤ k} \big)
= \text{det}\big( (M_{J(j),I(i)+1})_{1≤ i,j≤ k} \big),
```
and we obtain the components of ``ω̄^{α,J}=B_α φ^J`` in the basis
``\mathrm{d}x^I``
```math
ω̄_{I}^{α,J} = B_α \sum_{1≤l≤k+1} (-1)^{l+1} λ_{J(l)} \, m_I^{J\backslash l}.
```
The ``\binom{D}{k}\binom{N}{k}`` coefficients ``\{m_I^{J}\}_{I,J}`` are
constant in ``T`` and are pre-computed from ``M`` in
`_compute_PmΛ_basis_coefficients!` at the creation of `PmLambdaBasis`
and stored in its field `m`.

Finally, the pseudocode to evaluate our basis ``ω̄`` of ``ℙ_r^-Λ^k(T)`` at
``\boldsymbol{x}`` is
```julia
compute λ(x)
compute B(x) = { Bα(λ(x)) } for all |α|=r-1

for (F, bubble_functions) in get_bubbles(b)
    for (w, α, α_id, J, sub_J_ids) in bubble_functions

        ω̄_w = 0 # ω̄^{α,J}
        for (l, J_sub_Jl_id) in enumerate(sub_J_ids)
            λ_j = λ[J[l]]
            m_J_l = m[J_sub_Jl_id] # is a coordinate vector for all I
            ω̄_w += -(-1)^l * λ_j * m_J_l
        end

        Bα = B[α_id]
        ω̄[w] = Bα * ω̄_w
    end
end
```

#### Bubble functions ``\mathring{ℙ}_rΛ^k``

The ``ℙ`` type bubble basis polynomials associated to a face ``F⊆T`` defined by
[2, Th. 6.1-2] -- where the basis function Eq. (8.3) replace Eq. (8.1) -- are
```math
\mathring{ℙ}_rΛ^k(T,F) = \text{span}\big\{ ω^{α,J}=B_α Ψ^{α,J} \quad\big|\quad
\ α∈\mathcal{I}_{r}^D,\ \#J=k,\ ⟦α⟧∪J=F,\ α_i=0 \text{ if } i< \text{min}(F
\backslash J) \big\},
```
where ``Ψ^{α,J}`` [2, Eq. (8.3)] are defined by
```math
Ψ^{α,J} = \underset{j∈J}{⋀} Ψ^{α,F(α,J),j}
\quad\text{and}\quad
Ψ^{α,F,j} = \mathrm{d}λ^j - \frac{α_j}{|α|}\sum_{l∈F}\mathrm{d}λ^l,
```
where ``F(α,J)=⟦α⟧∪J``. `get_bubbles(b::PLambdaBasis)` provides the bubbles of
`b`, their bubble function indices are `(w, α, α_id, J)` only (do NOT modify
them).

Again, we need their components in the Cartesian basis
``\mathrm{d}x^I``:
```math
Ψ^{α,F,j} = M_{j,i+1}\mathrm{d}x^i - \frac{α_j}{|α|}\sum_{l∈F}M_{l,i+1}\mathrm{d}x^i
= \big(M_{j,i+1} - \frac{α_j}{|α|}\sum_{l∈F}M_{l,i+1}\big)\mathrm{d}x^i
```
so
```math
Ψ^{α,F,j} = ψ_{i}^{α,F,j} \mathrm{d}x^i
\quad\text{where}\quad
ψ_{i}^{α,F,j} = M_{j,i+1} - \frac{α_j}{|α|}\sum_{l∈F}M_{l,i+1}
```
and

```math
Ψ^{α,J} = ψ_I^{α,J} \mathrm{d}x^I
\quad\text{where}\quad
ψ_I^{α,J} = \text{det}\big( (ψ_{i}^{α,F,j})_{i∈I,\,j∈J} \big).
```

Finally, the ``\binom{D}{k}`` components of ``ω^{α,J}=B_α Ψ^{α,J}`` in the
basis ``\mathrm{d}x^I`` are
```math
ω_{I}^{α,J} = B_α\, ψ_I^{α,J},
```
where the ``\binom{D+r}{k+r}\binom{r+k}{k}\binom{D}{k}
=\mathrm{dim}(ℙ_rΛ^k(T^D))\times\# (\{\mathrm{d}x^I\}_I)`` coefficients
``ψ_I^{α,J}`` depend only on ``T`` and are pre-computed in
`_compute_PΛ_basis_form_coefficient!` at the construction of `PLambdaBasis` and
stored in its field `Ψ`.

The pseudocode to evaluate our basis ``ω`` of ``ℙ_rΛ^k(T)`` at
``\boldsymbol{x}`` is
```julia
compute λ(x)
compute B(x) = { Bα(λ(x)) } for all |α|=r

for (F, bubble_functions) in get_bubbles(b)
    for (w, α, α_id) in bubble_functions
        Bα = B[α_id]
        ω[w] = Bα * Ψ[w]
    end
end
```

#### Gradient and Hessian of the coefficient vectors

Let us derive the formula for the gradient and hessian of the basis forms
coefficient vectors ``\{ω̄_{I}^{α,J}\}_I`` and ``\{ω_{I}^{α,J}\}_I``. We will
express them in function of the scalar Bernstein polynomial derivatives already
implemented by `BernsteinBasisOnSimplex`. They are only supported for scalar or
`VectorValue`'d bases (vector calculus style).

##### Coefficient vector ``\{ω_{I}^{α,J}\}_I``

Recall ``ω_{I}^{α,J} = B_α\, ψ_I^{α,J}``. The derivatives are easy to
compute because only ``B_α`` depends on ``\boldsymbol{x}``, leading to
```math
∂_q\, ω_{I}^{α,J} = ψ_I^{α,J}\; ∂_q B_α,\qquad\text{or}\qquad ∇ω^{α,J} = ∇B_α ⊗ ψ^{α,J}\\

∂_t∂_q\, ω_{I}^{α,J} = ψ_I^{α,J}\; ∂_t∂_q B_α\qquad\text{or}\qquad ∇∇ω^{α,J} = ∇∇B_α ⊗ ψ^{α,J}.
```
where ``∇`` and ``∇∇`` are the standard gradient and hessian operators and
``ψ^{α,J}`` is seen as a length-``\binom{D}{k}`` vector with no variance
(not as a ``k``-form, which is an order ``k`` covariant tensor).

##### Coefficient vector ``\{ω̄_{I}^{α,J}\}_I``

Recall ``ω̄_{I}^{α,J} = B_α \sum_{1≤l≤k+1} (-1)^{l+1} λ_{J(l)} \,
m_I^{J\backslash l}``, the derivatives are not immediate to compute because
both ``B_α`` and ``λ_{J(l)}`` depend on ``\boldsymbol{x}``, let us first use ``B_α
λ_{J(l)} = \frac{α_{J(l)} + 1}{|α|+1}B_{α+e(J,l)}`` where ``e(J,l) =
\big(δ_i^{J_l}\big)_{1≤ i≤ N}`` to write the coefficients in Bernstein form as
follows
```math
ω̄_{I}^{α,J} = B_α \sum_{1≤l≤k+1} (-1)^{l+1} λ_{J(l)} \, m_I^{J\backslash l} =
\frac{1}{r}\sum_{1≤l≤k+1} (-1)^{l+1} (α_{J(l)} +1)\ B_{α+e(J,l)}\, m_I^{J\backslash l},
```
where ``|α|+1`` was replaced with ``r``, the polynomial degree of
``ω̄_{I}^{α,J}``. As a consequence, for any Cartesian coordinate indices
``1≤ p,q≤ D``, we get
```math
∂_q ω̄_{I}^{α,J} = \frac{1}{r}\sum_{1≤l≤k+1} (-1)^{l+1} (α_{J(l)} +1) \ ∂_q B_{α+e(J,l)}\, m_I^{J\backslash l},\\
∂_t∂_q ω̄_{I}^{α,J} = \frac{1}{r}\sum_{1≤l≤k+1} (-1)^{l+1} (α_{J(l)} +1) \ ∂_t∂_q B_{α+e(J,l)}\, m_I^{J\backslash l}.
```
In tensor form, this is
```math
\mathrm{D}ω^{α,J} = \frac{1}{r}\sum_{1≤l≤k+1} (-1)^{l+1} (α_{J(l)} +1)\ \mathrm{D}\!B_{α+e(J,l)} ⊗ m^{J\backslash l}\\
```
where ``\mathrm{D}`` is ``∇`` or ``∇∇``, the standard gradient and hessian operators,
and ``ω̄^{α,J}`` is again seen as a length-``\binom{D}{k}`` vector with no
variance.

#### Exterior derivative of the basis forms (to be implemented)

The exterior derivative of a ``k``-form ``ω=ω_{\tilde{I}}\,\mathrm{d}x^{\tilde{I}}`` is the
``k\!+\!1``-form
```math
\mathrm{d}ω = ∂_i ω_{\tilde{I}}\, \mathrm{d}x^i ∧ \mathrm{d}x^{\tilde{I}}
\qquad\text{ where }\qquad \#\tilde{I} = k
```
We need to express ``\mathrm{d}ω`` in the basis of ``k+1`` forms. Let ``I``
such that ``\#I=k\!+\!1`` with ``k<D`` (otherwise ``\mathrm{d}ω=0``). Because
the exterior product is alternating, the coefficients that contribute to
``(\mathrm{d}ω)_{I}`` are ``∂_i ω_{\tilde{I}}`` for which ``i=I_q`` and ``\tilde{I}
= I\backslash \{I_q\}`` with ``1≤ q≤ k+1``, so one can deduce
```math
(\mathrm{d}ω)_I = \underset{1≤ q≤ k+1}{\sum} (-1)^{q-1}\ ∂_{I(q)} ω_{I\backslash q}.
```

##### Polynomial forms ``\mathrm{d}\,ω̄^{α,J}``

For all ``|α|=r\!-\!1``, ``\,\#J=k\!+\!1`` and ``\#I = k\!+\!1`` (with ``k\!<\!D``):
```math
(\mathrm{d}\,ω̄^{α,J})_I = \frac{1}{r}\underset{1≤ l≤ k+1}{\sum} (-1)^{l+1}(α_{J(l)}+1)
\underset{1≤ q≤ k+1}{\sum} (-1)^{q-1}\ m_{I\backslash q}^{J\backslash l}\ ∂_{I(q)} B_{α+e(J,l)}\;
```

##### Polynomial forms ``\mathrm{d}\,ω^{α,J}``

For all ``|α|=r``, ``\,\#J=k`` and ``\#I = k\!+\!1`` (with ``k\!<\!D``):
```math
(\mathrm{d}\,ω^{α,J})_I =
\underset{1≤ q≤ k+1}{\sum} (-1)^{q-1}\ ψ_{I\backslash q}^{α,J}\ ∂_{I(q)} B_α.
```

#### Hodge operator of the basis forms

The Hodge operator of the canonical basis forms ``\mathrm{d}x^I`` in an
Euclidean (Riemannian) space is
```math
\star \mathrm{d}x^I = \mathrm{sgn}\left(I\!*\!\bar{I}\,\right)\mathrm{d}x^{\bar{I}}
```
where ``\bar{I}`` is the complement of ``I`` in ``1\!:\!D``, that is the only
combination such that the concatenated permutation ``I\!*\!\bar{I}`` is a
permutation of ``1\!:\!D``. ``\bar{I}`` is implemented by
[`_complement(I)`](@ref _complement). [`_combination_sign(I)`](@ref
_combination_sign) computes the sign of ``I\!*\!\bar{I}``.



#### Analytical formulas in the reference simplex

In the reference simplex ``\hat{T}``, the vertices and thus coefficients of
``M`` are known at compile time, so the coefficients ``m_I^J`` and
``ψ_I^{α,J}`` in ``\hat{T}``, denoted by ``\hat{m}_I^J`` and
``\hat{ψ}_I^{α,J}`` respectively, could be hard-coded at compile time in
`@generated` functions to avoid storing them in the basis and accessing them at
runtime. Let us derive the formulas for them.


##### Coefficients ``\hat{m}_I^J``

It was shown in the Barycentric coordinates section above that
``M_{j,i+1} = δ_{i+1,j} - δ_{1j}``. Let ``\#I=k`` and ``\#J=k``. We need
to compute the determinant of the matrix
```math
\hat{M}_{IJ}=(δ_{I(i)+1,\,J(j)}-δ_{1,J(j)})_{1≤ i,j≤ k}.
```
Let us define:
- ``s=δ_1^{J_1}``, that indicates if ``\hat{M}_{IJ}`` contains a column of ``-1``,
- ``p = \text{min } \{j\,|\, I_j+1 ∉J\}`` where ``\text{min}\,∅=0``, the index of the first row of ``\hat{M}_{IJ}`` containing no ``1``. ``p=0`` if and only if ``\hat{M}_{IJ}`` is the identity matrix,
- ``n =\# \{\ i\ |\ i>s,\, J_i-1∉I\}``, the number of columns of zeros of ``\hat{M}_{IJ}``.

Then it can be shown that ``k-n`` is the rank of ``\hat{M}_{IJ}``, and that
```math
\hat{m}_I^J = \mathrm{det}(\hat{M}_{IJ}) = (-1)^{p(I,J)}δ_0^{n(I,J)},
```
where the dependency of ``p,n`` on ``I,J`` is made explicit, so in ``\hat{T}``,
there is
```math
ω̄_{I}^{α,J} = B_α \sum_{1≤l≤k+1} (-1)^{l+1} λ_{J(l)} \, \hat{m}_I^{J\backslash l}.
```

##### Coefficients ``\hat{ψ}_I^{α,J}``

The expression of ``ψ_{i}^{α,F,j}`` in ``\hat{T}`` is
```math
\hat{ψ}_{i}^{α,F,j} = M_{j,i+1} - \frac{α_j}{|α|}\sum_{l∈F}δ_{i+1,l} - δ_{1l}
= M_{j,i+1} + \big( δ_{1,F_0}-\sum_{l∈F}δ_{i+1,l} \big)\frac{α_j}{|α|}
```
leading to
```math
\hat{ψ}_I^{α,J}  = \mathrm{det}\Big(\hat{M}_{IJ} + u\,v^{\intercal}\Big)
\quad\text{where}\quad
u^i = δ_{1,F(0)}-\sum_{l∈F}δ_{I(i)+1,\,l}, \qquad v^j = \frac{α_{J(j)}}{|α|}.
```
We can use the following matrix determinant lemma:
```math
\mathrm{det}(\hat{M}_{IJ} + uv^\intercal) =
\mathrm{det}(\hat{M}_{IJ}) + v^\intercal\mathrm{adj}(\hat{M}_{IJ})u.
```
The determinant ``\mathrm{det}(\hat{M}_{IJ})=\hat{m}_I^{J}`` was computed
above, but ``\mathrm{adj}(\hat{M}_{IJ})``, the transpose of the cofactor matrix
of ``\hat{M}_{IJ}``, is also needed. Let ``s=δ_1^{J_1}``, ``n`` and ``p`` be
defined as above, and additionally define
- ``q = \text{min } \{j\,|\,j>p,\ I_j+1 ∉J\}``, the index of the second row of ``\hat{M}_{IJ}`` containing no ``1`` (``q=0`` if there isn't any),
- ``m = \text{min } \{i\,|\,i>s,\ J_i-1 ∉I\}``, the index of the first column of ``\hat{M}_{IJ}`` containing only zeros (``m=0`` if there isn't any).

Then the following table gives the required information to apply the matrix
determinant lemma and formulas for ``\hat{ψ}_I^{α,J}``
```math
\begin{array}{|c|c|c|c|c|}
\hline
s  & n & \mathrm{rank}\hat{M}_{IJ} & \mathrm{adj}\hat{M}_{IJ} & \hat{ψ}_I^{α,J} \\
\hline
\hline
0   & 0 & k   & δ_{ij}                         & 1 + u \cdot v\\
\hline
0   & 1 & k-1 & (-1)^{m+p}δ_i^m δ^p_j          & (-1)^{m+p}v^m u^p \\
\hline
1   & 0 & k   & (-1)^p(δ_{J(i),\,I(j)+1}-δ_{p,J(j)})& (-1)^p(1-u^p|v|+\underset{1≤ l<p}{\sum}v^{l+1}u^l + \underset{p<l≤ k}{\sum}v^{l}u^l) \\
\hline\hspace{1mm}
1   & 1 & k-1 & (-1)^{m+p+q}δ_i^m(δ^q_j-δ^p_j) & (-1)^{m+p+q}v^m(u^q-u^p) \\
\hline
0/1 & \geq 2 & ≤ k-2 & 0 & 0 \\
\hline
\end{array}
```
In this table, ``m``, ``p`` and ``q`` depend on ``I`` and ``J``, ``u``
depends on ``F`` and ``I``, and ``v`` depends on ``α`` and ``J``.

### Low level docstrings

```@docs
_compute_cart_to_bary_matrix
_cart_to_bary
_de_Casteljau_nD!
_downwards_de_Casteljau_nD!
_combination_index
_sorted_combinations
_basis_forms_components
_combination_sign
_complement
```

## References

[1] [M.J. Lai & L.L. Schumaker, Spline Functions on Triangulations, Chapter 2 - Bernstein–Bézier Methods for Bivariate Polynomials, pp. 18 - 61.](https://doi.org/10.1017/CBO9780511721588.003)

[2] [D.N. Arnold, R.S. Falk & R. Winther, Geometric decompositions and local bases for spaces of finite element differential forms, Computer Methods in Applied Mechanics and Engineering](https://doi.org/10.1016/j.cma.2008.12.017)

[3] [D.N. Arnold and A. Logg, Periodic Table of the Finite Elements, SIAM News, vol. 47 no. 9, November 2014.](https://www-users.cse.umn.edu/~arnold/papers/periodic-table.pdf)
