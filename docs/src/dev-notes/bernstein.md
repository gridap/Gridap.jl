
# Bernstein basis algorithms

## Notations

The univariate [`Bernstein`](@ref) polynomials forming a basis of ``\mathbb{P}_K``
are defined by
```math
B^K_{n}(x) = \binom{K}{n} x^n (1-x)^{K-n}\qquad\text{ for } 0\leq n\leq K.
```

The ``D``-multivariate Bernstein polynomials of degree ``K`` are defined by
```math
B^{D,K}_\alpha(x) = \binom{K}{\alpha} λ^\alpha\qquad\text{for } |\alpha|=K
```
where
- ``\alpha`` belongs to ``\{ \alpha\in \}`` where ``N = D+1``
- ``\binom{K}{\alpha} = \frac{K!}{\alpha_1 !\alpha_2 !\dots\alpha_N!}`` and ``|\alpha|=\sum_{1\leq i\leq N} \alpha_i``

The superscript ``D`` and ``K`` can be omitted because they are always
determined by ``\alpha`` using ``D=\#(α)-1`` and ``K=|α|``.
The set ``\{B_\alpha\}_{\#(\alpha)=D+1,\,|\alpha|=K}`` is a basis of
``\mathbb{P}^D_K``, implemented by [`BernsteinBasisOnSimplex`](@ref).

## Pascal triangle, pyramid and simplex



## Standard de Casteljau algorithm


## Downward de Casteljau algorithm


## Evaluation of Bernstein polynomials and their derivatives

```math
B_\alpha = \sum_{1\leq i\leq N} λ_i B_{α-e_i} = \sum_{(i,β)\in sub(α)} λ_i B_β
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

## References

[1] [M.J. Lai & L.L. Schumaker, Spline Functions on Triangulations, Chapter 2 - Bernstein–Bézier Methods for Bivariate Polynomials, pp. 18 - 61.](https://doi.org/10.1017/CBO9780511721588.003)
