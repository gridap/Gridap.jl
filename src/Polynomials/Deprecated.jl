"""
    num_terms(a::PolynomialBasis)

!!! warning
    Deprecated in favor of length(a).
"""
function num_terms end
@deprecate num_terms(a::PolynomialBasis) length(a)

@deprecate MonomialBasis{D}(args...) where D MonomialBasis(Val(D), args...)

"""
    PCurlGradMonomialBasis{D}(T,order::Int) where D

!!! warning
    Deprecated in favor of `FEEC_poly_basis(Val(D),T,order+1,D-1,:P⁻,Monomial; rotate_90=(D==2))`.
"""
struct PCurlGradMonomialBasis{D}
  function PCurlGradMonomialBasis end # prevents any instantiation
end
@deprecate PCurlGradMonomialBasis{D}(::Type{T},order::Int) where {T,D} FEEC_poly_basis(Val(D),T,order+1,D-1,:P⁻,Monomial; rotate_90=(D==2))

"""
    QGradMonomialBasis{D}(T,order::Int) where D

!!! warning
    Deprecated in favor of `FEEC_poly_basis(Val(D),T,order+1,1,:Q⁻,Monomial)`.
"""
struct QGradMonomialBasis{D}
  function QGradMonomialBasis end
end
@deprecate QGradMonomialBasis{D}(::Type{T},order::Int) where {T,D} FEEC_poly_basis(Val(D),T,order+1,1,:Q⁻,Monomial)

"""
    QCurlGradMonomialBasis{D}(T,order::Int) where D

!!! warning
    Deprecated in favor of `FEEC_poly_basis(Val(D),T,order+1,D-1,:Q⁻,Monomial; rotate_90=(D==2))`.
"""
struct QCurlGradMonomialBasis{D}
  function QCurlGradMonomialBasis end
end
@deprecate QCurlGradMonomialBasis{D}(::Type{T},order::Int) where {T,D} FEEC_poly_basis(Val(D),T,order+1,D-1,:Q⁻,Monomial; rotate_90=(D==2))

struct NedelecPrebasisOnSimplex{D}
  function NedelecPrebasisOnSimplex end
end
@deprecate NedelecPrebasisOnSimplex{D}(order::Int) where D NedelecPolyBasisOnSimplex{D}(Monomial, Float64, order) false

"""
    JacobiPolynomialBasis{D}(args...) where D

!!! warning
    Deprecated in favor of LegendreBasis(Val(D), args...).
"""
struct JacobiPolynomialBasis{D}
  function JacobiPolynomialBasis end
end
@deprecate JacobiPolynomialBasis{D}(args...) where D LegendreBasis(Val(D), args...)
