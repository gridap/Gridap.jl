"""
    num_terms(a::PolynomialBasis)

!!! warning
    Deprecated in favor of length(a).
"""
function num_terms end

@deprecate num_terms(a::PolynomialBasis) length(a)

@deprecate MonomialBasis{D}(args...) where D MonomialBasis(Val(D), args...)

"""
    PGradMonomialBasis{D}(args...) where D

!!! warning
    Deprecated in favor of PGradBasis(Monomial, Val(D), args...).
"""
struct PGradMonomialBasis{D}
  function PGradMonomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate PGradMonomialBasis{D}(args...) where D PGradBasis(Monomial, Val(D), args...) false

"""
    PCurlGradMonomialBasis{D}(args...) where D

!!! warning
    Deprecated in favor of PCurlGradBasis(Monomial, Val(D), args...).
"""
struct PCurlGradMonomialBasis{D}
  function PCurlGradMonomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate PCurlGradMonomialBasis{D}(args...) where D PCurlGradBasis(Monomial, Val(D), args...)

"""
    QGradMonomialBasis{D}(args...) where D

!!! warning
    Deprecated in favor of QGradBasis(Monomial, Val(D), args...).
"""
struct QGradMonomialBasis{D}
  function QGradMonomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate QGradMonomialBasis{D}(args...) where D QGradBasis(Monomial, Val(D), args...)

"""
    QCurlGradMonomialBasis{D}(args...) where D

!!! warning
    Deprecated in favor of QCurlGradBasis(Monomial, Val(D), args...).
"""
struct QCurlGradMonomialBasis{D}
  function QCurlGradMonomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate QCurlGradMonomialBasis{D}(args...) where D QCurlGradBasis(Monomial, Val(D), args...)

struct NedelecPrebasisOnSimplex{D}
  function NedelecPrebasisOnSimplex()
    @unreachable
    new{0}()
  end
end
@deprecate NedelecPrebasisOnSimplex{D}(args...) where D PGradBasis(Monomial,Val(D), args...) false

"""
    JacobiPolynomialBasis{D}(args...) where D

!!! warning
    Deprecated in favor of LegendreBasis(Val(D), args...).
"""
struct JacobiPolynomialBasis{D}
  function JacobiPolynomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate JacobiPolynomialBasis{D}(args...) where D LegendreBasis(Val(D), args...)
