module GridapSymbolicsExt

using Symbolics
using Gridap
using Gridap.TensorValues
using Gridap.Polynomials

using Gridap.Polynomials: multinomial  # Combinatorics, via Gridap.Polynomials
using Gridap.Polynomials: _sup_str

import Gridap.TensorValues: exterior_derivative
import Gridap.TensorValues: codifferential
import Gridap.TensorValues: lie_derivative
import Gridap.TensorValues: symbolic_coordinates
import Gridap.Polynomials: print_forms

# ============================================================
# Symbolic scalar multiplication (Symbolics.Num is not <:Real)
# ============================================================

function Base.:*(s::Symbolics.Num, ω::DifferentialFormValue{K,D,T,L}) where {K,D,T,L}
  DifferentialFormValue{K,D,Symbolics.Num,L}(map(x -> s*x, ω.data))
end
Base.:*(ω::DifferentialFormValue, s::Symbolics.Num) = s * ω

# ============================================================
# Symbolic coordinate variables
# ============================================================

@variables x¹ x² x³ x⁴ x⁵ x⁶ x⁷ x⁸ x⁹
const coords = [x¹, x², x³, x⁴, x⁵, x⁶, x⁷, x⁸, x⁹]

symbolic_coordinates(D::Integer) = Tuple(coords[1:D])

# ============================================================
# Symbolic exterior derivative
# Dispatches when T = Symbolics.Num (symbolic coefficients).
# ============================================================

function _sym_exterior_derivative(ω::DifferentialFormValue{K,D,Symbolics.Num,L0},
                                   diff_vars::AbstractVector) where {K,D,L0}
  L   = binomial(D, K)
  ds  = [Symbolics.jacobian([ω.data[i]], diff_vars) for i in 1:L]
  dfs = [DifferentialFormValue{1,D}(Tuple(ds[i])) for i in 1:L]
  kbs = [DifferentialFormValue{K,D}(ntuple(i -> i == j ? 1.0 : 0.0, L)) for j in 1:L]
  res = [dfs[i] ∧ kbs[i] for i in 1:L]
  dat = sum([collect(r.data) for r in res])
  Lout = binomial(D, K+1)
  DifferentialFormValue{K+1,D,Symbolics.Num,Lout}(Tuple(dat))
end

# Differentiates w.r.t. global Cartesian coords x¹,...,x^D
exterior_derivative(ω::DifferentialFormValue{K,D,Symbolics.Num,L0}) where {K,D,L0} =
  _sym_exterior_derivative(ω, coords[1:D])

# Differentiates w.r.t. explicit variables — works for barycentric λs or any symbols
exterior_derivative(ω::DifferentialFormValue{K,D,Symbolics.Num,L0},
                    vars::NTuple{D,Symbolics.Num}) where {K,D,L0} =
  _sym_exterior_derivative(ω, collect(vars))

# ============================================================
# Symbolic codifferential  δ : Ωᴷ(D) → Ωᴷ⁻¹(D)
#   δ = (-1)^{D(K-1)+1} ⋆ d ⋆   (flat Euclidean metric, symbolic coefficients)
# ============================================================

function codifferential(ω::DifferentialFormValue{K,D,Symbolics.Num}) where {K,D}
  @assert K >= 1 "codifferential requires K ≥ 1"
  sgn = iseven(D*(K-1) + 1) ? 1 : -1
  sgn * exterior_derivative(hodge_star(ω))  |> hodge_star
end

# ============================================================
# Lie derivative  L_v : Ωᴷ(D) → Ωᴷ(D)
#   Cartan's magic formula:  L_v ω = d(ι_v ω) + ι_v(dω)
#   (symbolic coefficients required for d)
# ============================================================

function lie_derivative(v::VectorValue{D,Symbolics.Num},
                        ω::DifferentialFormValue{0,D,Symbolics.Num}) where D
  # L_v f = ι_v(df)  (directional derivative of a scalar)
  interior_product(v, exterior_derivative(ω))
end

function lie_derivative(v::VectorValue{D,Symbolics.Num},
                        ω::DifferentialFormValue{K,D,Symbolics.Num}) where {K,D}
  dω      = exterior_derivative(ω)               # (K+1)-form
  ι_v_ω   = interior_product(v, ω)               # (K-1)-form
  exterior_derivative(ι_v_ω) + interior_product(v, dω)
end

# ============================================================
# print_forms: symbolic display of the PΛ and P⁻Λ bases as ambient barycentric
# differential forms (dλ¹,…,dλ^{D+1})
# ============================================================

# A barycentric coordinate carries an upper index, so the symbols are λ¹,…,λᴺ
# rather than Symbolics.variables(:λ, 1:N), which names them λ₁,…,λ_N.
_barycentric_symbols(N) = [Symbolics.variable(Symbol("λ", _sup_str(i))) for i in 1:N]

const _λ_pow_parens = r"\((λ[⁰¹²³⁴⁵⁶⁷⁸⁹]+)\^(\d+)\)"
const _λ_pow_bare   = r"(λ[⁰¹²³⁴⁵⁶⁷⁸⁹]+)\^(\d+)"

# The basis, not the value, is what knows its components are coefficients on the
# ambient dλ coframe, so it is the basis that asks for the dλ labelling.
#
# Replaces λ¹^2 with (λ¹)²
function _show_barycentric(out::IO, ω)
  s = sprint((io, x) -> show(io, MIME("text/plain"), x), ω;
             context = (:coordinates => :barycentric,))
  for pat in (_λ_pow_parens, _λ_pow_bare)
    s = replace(s, pat => function (matched)
                            c = match(pat, matched)
                            "(" * c[1] * ")" * _sup_str(parse(Int, c[2]))
                          end)
  end
  print(out, s)
end

function print_forms(b::BarycentricPΛBasis{D}, out::IO=stdout) where D
  N  = D + 1
  r  = get_order(b)
  k  = b.k
  r ≥ 1 || throw(ArgumentError(
    "print_forms needs r ≥ 1: the PᵣΛᵏ basis forms of order 0 are the constant dxᴵ, use print_indices"))
  λs = _barycentric_symbols(N)

  # The direction 1-forms the basis builds in the physical dx¹,…,dxᴰ frame,
  # here on the ambient dλ¹,…,dλᴺ coframe instead: column j holds
  # φ^{α,F,j} = dλʲ - (αⱼ/r) Σ_{l∈F} dλˡ, with α weighted as the flavor asks.
  function ambient_φ(F, α, r, flavor)
    s, ρ = flavor === :BMM ? (map(αj -> Int(αj > 0), α), count(>(0), α)) : (α, r)
    φ = zeros(Float64, N, N)
    for j in 1:N
      φ[j,j] += 1
      for i in F
        φ[i,j] -= s[j]/ρ
      end
    end
    φ
  end

  println(out, "BarycentricPΛBasis{D=$D, r=$r, k=$k, $(b.flavor)}: dim = $(length(b))  (as barycentric differential forms)")
  for (F, bubble_functions) in get_bubbles(b)
    for (w, α, _, J) in bubble_functions
      mono = multinomial(α...) * prod(λs[i]^α[i] for i in 1:N)
      φ    = ambient_φ(F, α, r, b.flavor)  # ambient dλ¹,…,dλᴺ coefficients
      # The direction k-form is the wedge of the columns J of φ
      Ψ = foldl((ω, j) -> ω ∧ DifferentialFormValue{1,N}(Tuple(φ[:,j])),
                J; init=DifferentialFormValue{0,N}((1.0,)))
      print(out, "[", rpad(w,3), "]  F=", rpad(join(F,","),8), " J=", rpad(join(J,","),7), " α=", rpad(string(Tuple(α)),12), "  ")
      _show_barycentric(out, mono * Ψ)
      println(out)
    end
  end
end

function print_forms(b::BarycentricPmΛBasis{D}, out::IO=stdout) where D
  N  = D + 1
  r  = get_order(b)
  k  = b.k
  λs = _barycentric_symbols(N)

  dλ(j) = DifferentialFormValue{1,N}(ntuple(i -> i == j ? 1.0 : 0.0, N))

  # The Whitney k-form of a bubble, on the ambient dλ¹,…,dλᴺ coframe rather than
  # the physical dx¹,…,dxᴰ one the basis stores in its `m`:
  #   φ^J = Σ_l (-1)^{l+1} λ^{J(l)} dλ^{J∖J(l)}
  function ambient_φ(J)
    sum(enumerate(J)) do (l, Jl)
      dλ_JsubJl = foldl((ω, j) -> ω ∧ dλ(j), (j for j in J if j != Jl);
                        init=DifferentialFormValue{0,N}((1.0,)))
      (iseven(l) ? -λs[Jl] : λs[Jl]) * dλ_JsubJl
    end
  end

  println(out, "BarycentricPmΛBasis{D=$D, r=$r, k=$k, $(b.flavor)}: dim = $(length(b))  (as barycentric differential forms)")
  for (F, bubble_functions) in get_bubbles(b)
    for (w, α, _, J) in bubble_functions
      # :BMM scales the Whitney form by the bare monomial λ^α, :AFW by Bα
      coeff = b.flavor === :BMM ? 1 : multinomial(α...)
      mono  = coeff * prod(λs[i]^α[i] for i in 1:N)
      print(out, "[", rpad(w,3), "]  F=", rpad(join(F,","),8), " J=", rpad(join(J,","),7), " α=", rpad(string(Tuple(α)),12), "  ")
      _show_barycentric(out, mono * ambient_φ(J))
      println(out)
    end
  end
end

end # module
