"""
    LoopHelpers

Internal helper module implementing the matrices for Stam's algorithm.
All these helpers allocate new plain `Matrix`.

J. Stam, "Evaluation of Loop Subdivision Surfaces", 1998
https://research.cs.wisc.edu/graphics/Courses/559-f2001/Papers/559-Papers/Stam-Loop.pdf

$(public_names_in_md(@__MODULE__))
"""
module LoopHelpers

using LinearAlgebra: Diagonal, diag, pinv, det, I, inv
using Gridap.Helpers: public_names_in_md

# Around an extraordinary vertex (EV) of valence `N`, vertices are grouped as
# `[EV, v2=ring_1, …, v(N+1)=ring_N, outer_1, …, outer_5]` (`K = N+6` total), the same
# convention used by `LoopPatchVerticesMap`/Stam's Fig. 2: the EV, its `N`-cycle
# of ring vertices, then 5 more vertices completing the patch. One further
# subdivision level needs 6 more vertices, for `M = N+12` total (Stam's Fig. 3).

# This is the graph of the extended neighborhood around the EV 1, after applying the extended subdivision matrix A_bar
"""

```plain
5------6 ···· N------N+6----N+12
|\\     |     /|     /|     /|
| \\    |   /  |   /  |   /  |
|   \\  |  /   |  /   |  /   |
|     \\|/     |/     |/     |
4------1------N+1----N+5----N+11
|     /|     /|     /|     /
|   /  |   /  |   /  |   /
|  /   |  /   |  /   |  /
|/     |/     |/     |/
3------2------N+2----N+10
|     /|     /|     /
|   /  |   /  |   /
|  /   |  /   |  /
|/     |/     |/
N+4----N+3----N+7
|     /|     /
|   /  |   /
|  /   |  /
|/     |/
N+9----N+8
```
"""

"""
    α(N) = 5/8 - (3 + 2*cos(2π/N))^2/64

Loop's extraordinary-vertex weight for valence `N`.
"""
α(N) = 5/8 - (3 + 2*cos(2π/N))^2/64

"""
    f(k, N) = 3/8 + (2/8)*cos(2π*k/N)

Eigenvalue of the circulant ring part of `S(N)` at ring frequency `k = 1,…,N-1`.
"""
f(k, N) = 3/8 + (2/8)*cos(2π*k/N)

"""
    Φ(N, n, k) -->  Λᴺ * A_tildeᵀ = Λᴺ * (Pₖ * A_bar * V)ᵀ
"""
Φ(N, n, k) = Λn(N,n) * transpose(A_tilde(N,k))

"""
    A_tilde(N, k) -->  Pₖ * A_bar * V
"""
A_tilde(N, k) = P(N,k) * A_bar(N) * V(N)

"""
    A(N)

`K×K` (`K=N+6`) subdivision matrix
```plain
A = [ S     0  ]
    [ S11  S12 ]
```
"""
function A(N)
  K = N+6
  M = zeros(K, K)
  M[1:N+1, 1:N+1] = S(N)
  M[N+2:K, 1:N+1] = S11(N)
  M[N+2:K, N+2:K] = S12()
  M
end

"""
    A_bar(N)

`(K+6)×K` extended subdivision matrix
```plain
A = [  S    0  ]
    [ S11  S12 ]
    [ S12  S13 ]
```
"""
function A_bar(N)
  K, Mtot = N+6, N+12
  M = zeros(Mtot, K)
  M[1:K, :]          = A(N)
  M[K+1:Mtot, 1:N+1] = S21(N)
  M[K+1:Mtot, N+2:K] = S22()
  M
end

"""
    P(N, k)

`12×(N+12)` picking matrix selecting the 12 control vertices of the
`k`-th (`k=1,2,3`) regular sub-patch

For `N=3`, `P(N,2)` selects the same vertex at two different Fig. 1 positions
(`v3` and `vN` coincide) — a genuine degeneracy of the minimum valence.
"""
function P(N, k)
  cols = if k == 1
    (N+9, N+8, N+4, N+3, N+7, 3, 2, N+2, N+10, 1, N+1, N+5)
  elseif k == 2
    (N+6, N, N+5, N+1, 1, N+10, N+2, 2, 3, N+7, N+3, N+4)
  elseif k == 3
    (N+3, N+7, 2, N+2, N+10, 1, N+1, N+5, N+11, N, N+6, N+12)
  else
    throw(ArgumentError("k must be 1, 2 or 3"))
  end
  M = zeros(12, N+12)
  for (row,col) in enumerate(cols)
    M[row,col] = 1.0
  end
  M
end

"""
    S(N)

`(N+1)×(N+1)` subdivision matrix of the EV's 1-ring neighborhood.
"""
function S(N)
  a, b, c, d = 1 - α(N), α(N)/N, 3/8, 1/8
  M = zeros(N+1, N+1)
  M[1,1] = a
  M[1,2:end] .= b
  for i in 1:N
    r = 1+i
    M[r,1]               += c
    M[r,1+i]             += c
    M[r,1+mod1(i-1,N)]   += d
    M[r,1+mod1(i+1,N)]   += d
  end
  M
end

"""
    S_hat(N)

Fourier transform of [`S(N)`](@ref)
"""
function S_hat(N)
  a, b, c, d = 1 - α(N), α(N)/N, 3/8, 1/8
  M = zeros(N+1, N+1)
  M[1,1], M[1,2] = a, N*b
  M[2,1], M[2,2] = c, c + 2d
  for k in 1:N-1
    M[2+k,2+k] = f(k,N)
  end
  M
end

"""
    Σ(N)

Diagonal matrix of the eigenvalues of `S(N)`, ordered to match the columns of `U0(N)`.
"""
function Σ(N)
  λ = zeros(N+1)
  λ[1], λ[2] = 1.0, 5/8 - α(N)
  col = 3
  for k in 1:(N-1)÷2
    λ[col], λ[col+1] = f(k,N), f(k,N)
    col += 2
  end
  iseven(N) && (λ[col] = f(N÷2, N))
  Diagonal(λ)
end

"""
    U0(N)

Real eigenvectors of `S(N)` (as columns), ordered to match `Σ(N)`.
"""
function U0(N)
  M = zeros(N+1, N+1)
  M[:,1] .= 1.0
  a, b, λ2 = 1 - α(N), α(N)/N, 5/8 - α(N)
  M[1,2] = -N*b/(a - λ2)
  M[2:end,2] .= 1.0
  col = 3
  for k in 1:(N-1)÷2
    for i in 1:N
      M[1+i,col]   = cos(2π*k*(i-1)/N)
      M[1+i,col+1] = sin(2π*k*(i-1)/N)
    end
    col += 2
  end
  if iseven(N)
    k = N÷2
    for i in 1:N
      M[1+i,col] = cos(2π*k*(i-1)/N)
    end
  end
  M
end

"""
    S11(N)

`5×(N+1)` block coupling the EV 1-ring to the 5 subdivided outer vertices of the `K=N+6` patch.
"""
function S11(N)
  M = zeros(5, N+1)
  M[:,1]   .+= [2, 1, 2, 1, 2]  # EV
  M[:,2]   .+= [6, 10, 6, 1, 0] # v2
  M[:,3]   .+= [0, 1, 6, 0, 0]  # v3
  M[:,N]   .+= [0, 0, 0, 1, 6]  # vN
  M[:,N+1] .+= [6, 1, 0, 10, 6] # vN+1
  M ./ 16
end

"""
    S12()

Fixed `5×5` block coupling the 5 subdivided outer vertices of the `K=N+6` patch to themselves.
"""
S12() = [2 0 0 0 0
         1 1 1 0 0
         0 0 2 0 0
         1 0 0 1 1
         0 0 0 0 2] ./ 16

"""
    S21(N)

`6×(N+1)` block coupling the EV 1-ring to the 6 subdivided 2-ring vertices of the K+6 patch.
"""
function S21(N)
  M = zeros(6, N+1)
  M[:,2]   .+= [3, 3, 3, 1, 0, 0] # # v2
  M[:,3]   .+= [0, 0, 1, 0, 0, 0] # # v3
  M[:,N]   .+= [0, 0, 0, 0, 0, 1] # # vN
  M[:,N+1] .+= [1, 0, 0, 3, 3, 3] # # vN+1
  M ./ 8
end

"""
    S22()

Fixed `6×5` block coupling the EV 1-ring to the 6 subdivided 2-ring Vertices.
Fixed `6×5` block coupling the 5 outer vertices
of the `K=N+6` patch to the 6 second-level vertices.
"""
S22() = [3 1 0 0 0
         1 3 1 0 0
         0 1 3 0 0
         3 0 0 1 0
         1 0 0 3 1
         0 0 0 1 3] ./ 8

"""
    Δ()

Diagonal matrix of the eigenvalues of `S12()`, ordered to match the columns of `W1()`.
"""
Δ() = Diagonal([1/8, 1/8, 1/8, 1/16, 1/16])

"""
    W1()

Eigenvectors of `S12()` (as columns), ordered to match `Δ()`.
"""
W1() = Float64[0 -1 1 0 0
               1 -1 1 0 1
               1  0 0 0 0
               0  0 1 1 0
               0  1 0 0 0]

"""
    U1(N)

Extension of `U0(N)`'s eigenvectors to the 5 outer vertices, solving
`U1*Σ(N) - S12()*U1 = S11(N)*U0(N)` column by column (Stam's Eq. 4).

Two kinds of eigenvalues of `Σ(N)` can coincide with an eigenvalue of `S12()`,
making that column's `5×5` system singular:
- `N` even, last column (`Σ(N)[N+1,N+1] = f(N/2,N) = 1/8`, multiplicity 3 in
  `S12()`): Stam gives the closed-form fix `(0,8,0,-8,0)`.
- Any other coincidence (e.g. `N=3`, where `Σ(N)[2,2] = 1/16`) is a genuine
  null-space ambiguity Stam's paper does not resolve; the minimum-norm
  solution is used instead.
"""
function U1(N)
  Σv, s12, rhs = diag(Σ(N)), S12(), S11(N)*U0(N)
  W = zeros(5, N+1)
  for j in 1:N+1
    if iseven(N) && j == N+1
      W[:,j] = [0.0, 8.0, 0.0, -8.0, 0.0]
      continue
    end
    Mj = Σv[j]*I - s12
    W[:,j] = abs(det(Mj)) < 1e-10 ? pinv(Mj)*rhs[:,j] : Mj\rhs[:,j]
  end
  W
end

"""
    V(N)

Eigenvector matrix of `A(N)` (as columns):
```plain
V = [ U0   0  ]
    [ U1  W1  ]
```
For `N=3`, column 2 (`U0(3)[:,2]` on top) is instead paired with a
generalized eigenvector (see [`_jordan_fix_N3`](@ref)), and columns
`N+1+4`/`N+1+5` (`Δ`'s two `1/16` entries) hold the corresponding genuine
eigenvector and its complement, rather than the raw `W1` columns.
"""
function V(N)
  K = N+6
  M = zeros(K, K)
  M[1:N+1, 1:N+1] = U0(N)
  M[N+2:K, 1:N+1] = U1(N)
  M[N+2:K, N+2:K] = W1()
  if N == 3
    u1, wa, wb = _jordan_fix_N3()
    M[N+2:K, 2]     = u1
    M[N+2:K, N+1+4] = wa
    M[N+2:K, N+1+5] = wb
  end
  M
end

# `A(3)` is analytically defective: the shared eigenvalue `1/16` (`σ(3)[2]`,
# coinciding with two of `δ()`'s entries) has algebraic multiplicity 3 but
# geometric multiplicity 2, so `eq. 4` has no solution for `u1(3)`'s column 2.
# `v`/`λ` below implements the generlized eigenvectors for the jordan block of 1/16.
_jordan_fix_N3() = ([-9, 0, -9, 0, -9], [0, -15/32, 0, -15/32, 0], [0, -1, 0, 1, 0])

"""
    V_inv(N)

```plain
V⁻¹ = [ U0⁻¹            0   ]  where V = [ U0   0  ]
      [ -W1⁻¹*U1*U0⁻¹  W1⁻¹ ]            [ U1  W1  ]
```
"""
function V_inv(N)
  K = N+6
  V_N = V(N)
  U0_N = V_N[1:N+1, 1:N+1]
  U1_N = V_N[N+2:K, 1:N+1]
  W1_N = V_N[N+2:K, N+2:K]
  U0i, W1i = inv(U0_N), inv(W1_N)
  M = zeros(K, K)
  M[1:N+1, 1:N+1] = U0i
  M[N+2:K, 1:N+1] = -W1i*U1_N*U0i
  M[N+2:K, N+2:K] = W1i
  M
end

"""
    Λ(N)

(Generalized) eigenvalue matrix of `A(N)`: `A = V*Λ*V_inv`.

Equal to `blockdiag(Σ(N), Δ())` for all `N` except `N=3`, where `Λ` has a jordan
block with one extra off-diagonal `1` at `Λ(3)[N+1+4, 2]`.
"""
function Λ(N)
  K = N+6
  M = zeros(K, K)
  M[1:N+1, 1:N+1] = Σ(N)
  M[N+2:K, N+2:K] = Δ()
  N == 3 && (M[N+1+4, 2] = 1.0)
  M
end

"""
    L(i, N)

`L(i, V_inv(N))`.
"""
L(i, N::Integer) = L(i, V_inv(N))

"""
    L(i, V_inv)

`Lᵢ = V_inv[ordᵢ, :]`, `ord = sortperm(eigenvalues, rev=true)`: `i=1` is the
limit position mask, `i=2,3` are the tangent-plane masks.
"""
function L(i, Vi::AbstractMatrix)
  N = size(Vi, 1) - 6
  λ = vcat(Σ(N).diag, Δ().diag)
  ord = sortperm(λ, rev=true)
  Vi[ord[i], :]
end

"""
    Λn(N, n) --> Λ(N)^n
"""
function Λn(N, n)
  K = N+6
  d = vcat(Σ(N).diag, Δ().diag)
  M = zeros(K, K)
  for i in 1:K
    M[i,i] = d[i]^n
  end
  if N == 3
    λ = Δ()[4,4]
    M[N+1+4, 2] = n == 0 ? 0.0 : n*λ^(n-1)
  end
  M
end


export α, f, S, S_hat, Σ, U0, S11, S12, S21, S22, Δ, W1, U1, A, A_bar, P, V,
       V_inv, Λ, Λn, A_tilde, Φ, L

end # module LoopHelpers
