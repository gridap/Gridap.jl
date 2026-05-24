using Printf
using Gridap
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Polynomials: BernsteinBasisOnSimplex, BarycentricPΛBasis,
                          BarycentricPmΛBasis, CartProdPolyBasis,
                          CompWiseTensorPolyBasis, ModalC0
using Gridap.Fields: evaluate
using LinearAlgebra: norm, dot

# ─────────────────────────────────────────────────────────────────────────────
# Core function
# ─────────────────────────────────────────────────────────────────────────────

"""
    evaluate_at_permutations(poly::Polytope, basis; quad_degree=nothing)
        -> (vertex_perms, E_ref, E_perms)

Evaluates `basis` at a fixed set of quadrature points on `poly`, then at the
corresponding permuted points for each admissible vertex permutation of `poly`.

Returns:
- `vertex_perms`: the list from `get_vertex_permutations(poly)`
- `E_ref`:   `(nq × nbasis)` matrix, `E_ref[q, i]  = φᵢ(xq[q])`
- `E_perms`: one `(nq × nbasis)` matrix per permutation σ,
             `E_perms[k][q, i] = φᵢ(Mσ(xq[q]))` where Mσ is the affine
             self-map of `poly` induced by σ.

If a basis function is equivariant, then column `i` of `E_perms[k]` equals
column `π(i)` of `E_ref` (up to a scalar for signed/scaled cases).
"""
function evaluate_at_permutations(poly::Polytope, basis; quad_degree=nothing)
  deg = isnothing(quad_degree) ? 2*get_order(basis) + 2 : quad_degree
  quad = Quadrature(poly, deg)
  xq   = get_coordinates(quad)

  E_ref = evaluate(basis, xq)

  lin_sf = get_shapefuns(LagrangianRefFE(Float64, poly, 1))
  verts  = get_vertex_coordinates(poly)

  vertex_perms = get_vertex_permutations(poly)

  E_perms = map(vertex_perms) do σ
    pverts = verts[invperm(σ)]
    xq_σ   = evaluate(linear_combination(pverts, lin_sf), xq)
    evaluate(basis, xq_σ)
  end

  vertex_perms, E_ref, E_perms
end

# ─────────────────────────────────────────────────────────────────────────────
# Analysis helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    analyse_column(col_σ, E_ref; tol=1e-10)

For a single column `col_σ` of a permuted evaluation matrix, find which column
`j` of `E_ref` it is most proportional to, the proportionality constant `c`,
and the residual `‖col_σ − c·E_ref[:,j]‖`.

Returns `(j, c, residual)`.
"""
# Works for both scalar (Float64) and vector-valued (VectorValue) columns.
_col_dot(a, b) = sum(ai ⊙ bi for (ai, bi) in zip(a, b))
_col_nrm(a)   = sqrt(real(_col_dot(a, a)))

# Flatten a column to a plain Float64 vector (needed for the 2-col least squares).
function _flatten(col)
  T = eltype(col)
  T <: Real && return Float64.(col)
  # VectorValue or similar MultiValue: unpack each entry's components
  Float64[xi for v in col for xi in Tuple(v)]
end

function analyse_column(col_σ, E_ref; tol=1e-10)
  n = size(E_ref, 2)
  best_j, best_c, best_res = 0, NaN, Inf
  for j in 1:n
    col_j = E_ref[:, j]
    nrm   = _col_nrm(col_j)
    if nrm < tol
      if _col_nrm(col_σ) < tol
        best_j, best_c, best_res = j, 1.0, 0.0
        break
      end
      continue
    end
    c   = _col_dot(col_σ, col_j) / nrm^2
    res = _col_nrm(col_σ .- c .* col_j)
    if res < best_res
      best_j, best_c, best_res = j, real(c), res
    end
  end
  best_j, best_c, best_res
end

"""
    analyse_column_2(col_σ, E_ref; tol=1e-10)

Like `analyse_column` but searches over all pairs of columns in `E_ref`,
solving the 2-column least-squares problem

    min_{c1,c2} ‖col_σ − c1·col_ref[:,j1] − c2·col_ref[:,j2]‖

Returns `(j1, j2, c1, c2, residual)` for the best-fitting pair.
Works for both scalar and vector-valued (VectorValue) entries by flattening
to plain Float64 vectors before solving.
"""
function analyse_column_2(col_σ, E_ref; tol=1e-10)
  n  = size(E_ref, 2)
  b  = _flatten(col_σ)
  cols = [_flatten(E_ref[:, j]) for j in 1:n]

  best_j1, best_j2, best_c1, best_c2, best_res = 0, 0, NaN, NaN, Inf
  for j1 in 1:n, j2 in j1+1:n
    A   = hcat(cols[j1], cols[j2])
    x   = A \ b          # least-squares via QR
    res = norm(b .- A * x)
    if res < best_res
      best_j1, best_j2, best_c1, best_c2, best_res = j1, j2, x[1], x[2], res
    end
  end
  best_j1, best_j2, best_c1, best_c2, best_res
end

"""
    classify_permutation(E_ref, E_σ; tol=1e-8) -> Symbol

Classify the relationship between `E_ref` and `E_σ`:
- `:pure`   — each column of E_σ equals a column of E_ref (pure permutation)
- `:signed` — each column equals ±1 times a column of E_ref (signed permutation)
- `:scaled` — each column is proportional to some column of E_ref (scaled permutation)
- `:none`   — some column is not proportional to any column of E_ref
"""
function classify_permutation(E_ref, E_σ; tol=1e-8)
  n = size(E_ref, 2)
  all_pure, all_signed = true, true
  scale = _col_nrm(E_ref[:, 1])
  for i in 1:n
    _, c, res = analyse_column(E_σ[:, i], E_ref; tol)
    res > tol * (scale + tol)  && return :none
    abs(abs(c) - 1.0) > tol   && (all_pure = false; all_signed = false)
    abs(c - 1.0) > tol        && (all_pure = false)  # c = -1 is signed, not pure
  end
  all_pure ? :pure : all_signed ? :signed : :scaled
end

# ─────────────────────────────────────────────────────────────────────────────
# Pretty-printing
# ─────────────────────────────────────────────────────────────────────────────

"""
    print_permutation_analysis(poly, basis, label=""; tol=1e-8, show_matrices=false)

Run `evaluate_at_permutations` and pretty-print a table showing, for each
vertex permutation σ, how each basis function transforms:

    φᵢ(Mσ(x)) ≈ c · φⱼ(x)

If `show_matrices=true`, also prints the full `E_ref` and each `E_σ` matrix
(rows = quadrature points, columns = basis functions) so you can inspect the
raw evaluations directly.
"""
function print_permutation_analysis(poly::Polytope, basis, label="";
                                    tol=1e-8, show_matrices=false)
  vertex_perms, E_ref, E_perms = evaluate_at_permutations(poly, basis)
  n = size(E_ref, 2)

  lbl = isempty(label) ? "" : " [$label]"
  println("━"^72)
  println("  Polytope: $poly on $(typeof(basis).name.name)")
  println("  n_basis=$(n), n_perms=$(length(vertex_perms))$lbl")
  println("━"^72)

  if show_matrices
    println("  E_ref  (rows=quad pts, cols=basis funs):")
    _print_matrix(E_ref)
    println()
  end

  for (σ, E_σ) in zip(vertex_perms, E_perms)
    kind = classify_permutation(E_ref, E_σ; tol)
    kind_str = Dict(:pure=>"\e[32mpure perm\e[0m", :signed=>"\e[33msigned perm\e[0m",
                    :scaled=>"\e[35mscaled perm\e[0m", :none=>"\e[31mNOT a perm\e[0m")[kind]
    println("  σ = $σ  →  $kind_str")

    if kind === :pure
      perm = [analyse_column(E_σ[:, i], E_ref; tol)[1] for i in 1:n]
      println("       π = $perm")
    elseif kind === :signed
      for i in 1:n
        j, c, _ = analyse_column(E_σ[:, i], E_ref; tol)
        sgn = c > 0 ? "+" : "-"
        println("     φ[$i](Mσ·) = $sgn φ[$j]")
      end
    else
      for i in 1:n
        j, c, res = analyse_column(E_σ[:, i], E_ref; tol)
        scale = _col_nrm(E_ref[:, 1])
        bad   = res > tol * (scale + tol)

        if !bad
          # Single-column fit is fine
          c_str = @sprintf("%.4f", c)
          r_str = res < 1e-12 ? "≈0" : @sprintf("%.2e", res)
          println("     φ[$i](Mσ·) ≈ $c_str · φ[$j]   (res=$r_str)")
        else
          # Try 2-column combination
          j1, j2, c1, c2, res2 = analyse_column_2(E_σ[:, i], E_ref; tol)
          good2 = n > 1 && res2 < tol * (scale + tol)
          if good2
            c1s = @sprintf("%.4f", c1)
            c2s = @sprintf("%.4f", c2)
            r2s = res2 < 1e-12 ? "≈0" : @sprintf("%.2e", res2)
            println("     φ[$i](Mσ·) ≈ $c1s·φ[$j1] + $c2s·φ[$j2]   (res=$r2s)")
          else
            # Neither fit worked — show best 1-col with flag
            c_str = @sprintf("%.4f", c)
            r_str = @sprintf("%.2e", res)
            println("  !  φ[$i](Mσ·) ≈ $c_str · φ[$j]   (1-col res=$r_str, 2-col res=$(@sprintf("%.2e",res2)))")
          end
        end
      end
    end

    if show_matrices
      println("  E_σ  (σ=$σ):")
      _print_matrix(E_σ)
    end
    println()
  end
end

# Format a single evaluation value for matrix display.
_fmt(x::Real)                   = @sprintf("%8.4f", x)
_fmt(x::VectorValue)            = "[" * join(_fmt.(Tuple(x)), ", ") * "]"
_fmt(x)                         = string(x)

function _print_matrix(E::AbstractMatrix; indent="    ")
  nq, nb = size(E)
  for q in 1:nq
    row = join((_fmt(E[q, i]) for i in 1:nb), "  ")
    println(indent, "q$q: ", row)
  end
end

# ─────────────────────────────────────────────────────────────────────────────
# Examples
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "═"^72)
println("  BASIS PERMUTATION ANALYSIS")
println("  Confirming pure / signed / scaled / none cases")
println("═"^72 * "\n")

# ── 1. BernsteinBasisOnSimplex on SEGMENT (pure perm) ─────────────────────
print_permutation_analysis(SEGMENT,
  BernsteinBasisOnSimplex{1}(Float64, 2),
  "BernsteinBasisOnSimplex{1} order 2  [expect: pure]")

# ── 2. BarycentricPΛBasis k=0 on SEGMENT (pure perm) ─────────────────────
print_permutation_analysis(SEGMENT,
  BarycentricPΛBasis{1}(Float64, 2, 0),
  "BarycentricPΛBasis{1} k=0 order 2  [expect: pure]")

# ── 3. BarycentricPmΛBasis k=0 on SEGMENT, order 2 (pure perm up to r=2) ─
print_permutation_analysis(SEGMENT,
  BarycentricPmΛBasis{1}(Float64, 2, 0),
  "BarycentricPmΛBasis{1} k=0 order 2  [expect: pure for r≤2]")

# ── 4. BarycentricPmΛBasis k=0 on SEGMENT, order 3 (scaled, c=2) ─────────
print_permutation_analysis(SEGMENT,
  BarycentricPmΛBasis{1}(Float64, 3, 0),
  "BarycentricPmΛBasis{1} k=0 order 3  [expect: scaled c=2]")

# ── 5. BarycentricPΛBasis k=1 on SEGMENT  ─────────────────────────────────
#    (RT TRI cell moment  /  ND1 TRI cell moment)
print_permutation_analysis(SEGMENT,
  BarycentricPΛBasis{1}(Float64, 1, 1),
  "BarycentricPΛBasis{1} k=1 order 1  [expect: signed]")

# ── 6. BarycentricPΛBasis k=0 on TRI (RT TET facet, ND1 TET edge basis) ──
print_permutation_analysis(TRI,
  BarycentricPΛBasis{2}(Float64, 2, 0),
  "BarycentricPΛBasis{2} k=0 order 2  [expect: pure]")

# ── 7. BarycentricPΛBasis k=1 on TRI (ND1 TET triangular face moments) ───
print_permutation_analysis(TRI,
  BarycentricPΛBasis{2}(Float64, 1, 1),
  "BarycentricPΛBasis{2} k=1 order 1  [expect: signed]")

# ── 8. CartProdPolyBasis{1,ModalC0} order 2 on SEGMENT (RT/ND QUAD edge) ──
print_permutation_analysis(SEGMENT,
  CartProdPolyBasis(ModalC0, Val(1), Float64, 2),
  "CartProdPolyBasis{1,ModalC0} order 2  [expect: signed]")

# ── 9. CartProdPolyBasis{1,Bernstein} order 2 on SEGMENT (pure perm) ──────
print_permutation_analysis(SEGMENT,
  CartProdPolyBasis(Bernstein, Val(1), Float64, 2),
  "CartProdPolyBasis{1,Bernstein} order 2  [expect: pure]")

# ── 10. CartProdPolyBasis{2,Bernstein} on QUAD (pure perm) ────────────────
print_permutation_analysis(QUAD,
  CartProdPolyBasis(Bernstein, Val(2), Float64, 2),
  "CartProdPolyBasis{2,Bernstein} Q₂ on QUAD  [expect: pure]")
