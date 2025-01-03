@deprecate num_terms(a::PolynomialBasis) length(a)

@deprecate MonomialBasis{D}(args...) where D MonomialBasis(Val(D), args...)

struct PGradMonomialBasis{D}
  function PGradMonomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate PGradMonomialBasis{D}(args...) where D PGradBasis(Monomial, Val(D), args...) false

struct PCurlGradMonomialBasis{D}
  function PCurlGradMonomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate PCurlGradMonomialBasis{D}(args...) where D PCurlGradBasis(Monomial, Val(D), args...)

struct QGradMonomialBasis{D}
  function QGradMonomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate QGradMonomialBasis{D}(args...) where D QGradBasis(Monomial, Val(D), args...)

struct QCurlGradMonomialBasis{D}
  function QCurlGradMonomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate QCurlGradMonomialBasis{D}(args...) where D QCurlGradBasis(Monomial, Val(D), args...)

struct NedelecPreBasisOnSimplex{D}
  function NedelecPreBasisOnSimplex()
    @unreachable
    new{0}()
  end
end
@deprecate NedelecPreBasisOnSimplex{D}(args...) where D NedelecPolyBasisOnSimplex{D}(args...) false

struct JacobiPolynomialBasis{D}
  function JacobiPolynomialBasis()
    @unreachable
    new{0}()
  end
end
@deprecate JacobiPolynomialBasis{D}(args...) where D LegendreBasis(Val(D), args...)
