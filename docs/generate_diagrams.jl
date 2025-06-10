using Kroki

macro export_kroki(name, str)
  return quote
    write("src/assets/"*string($name)*".svg", render($str, "svg"))
  end
end

# Plantuml uml diagrams syntax: https://plantuml.com/class-diagram#49b7759afaffc066

@export_kroki :poly_1 plantuml"""
  @startuml
  skinparam groupInheritance 2
  hide empty members

  together {
    abstract Field
    abstract Polynomial {
      +isHierarchical
    }
  }
  Field <|-left- Polynomial

  together {
    struct Monomial
    struct Legendre
    struct Chebyshev
    struct ModalC0
    struct Bernstein
  }

  Polynomial <|-- Monomial
  Polynomial <|-- Legendre
  Polynomial <|-- Chebyshev
  Polynomial <|-- ModalC0
  Polynomial <|-- Bernstein

  @enduml
  """

@export_kroki :poly_2 plantuml"""
  @startuml
  skinparam groupInheritance 2
  hide empty members

  together {
    abstract "AbstractArray{<:Polynomial}" as a1
    abstract PolynomialBasis {
      +get_order
      +return_type
    }
  }
  a1 <|-left- PolynomialBasis

  together {
    struct CartProdPolyBasis {
      +get_exponents
      +get_orders
    }
    struct CompWiseTensorPolyBasis
    struct RaviartThomasPolyBasis
    struct NedelecPolyBasisOnSimplex
    struct BernsteinBasisOnSimplex
    struct ModalC0Basis {
      +get_orders
    }
  }

  PolynomialBasis <|-- CartProdPolyBasis
  PolynomialBasis <|-- CompWiseTensorPolyBasis
  PolynomialBasis <|-- RaviartThomasPolyBasis
  PolynomialBasis <|-- NedelecPolyBasisOnSimplex
  PolynomialBasis <|-- BernsteinBasisOnSimplex
  PolynomialBasis <|-- ModalC0Basis

  object "(<:Polynomial)Basis" as m1
  object "QGrad[<:Polynomial]Basis\nQCurlGrad[<:Polynomial]Basis" as m2
  object "PCurlGrad[<:Polynomial]Basis" as m4
  object "PGradMonomialBasis" as m5
  CartProdPolyBasis <-down- m1
  CompWiseTensorPolyBasis <-down- m2
  RaviartThomasPolyBasis <-down- m4
  NedelecPolyBasisOnSimplex <-down- m5

  @enduml
  """

