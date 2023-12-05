###############
# TableauType #
###############
"""
    abstract type TableauType <: GridapType end

Trait that indicates whether a tableau is explicit, implicit or
implicit-explicit.
"""
abstract type TableauType <: GridapType end

"""
    struct ExplicitTableau <: TableauType end

Tableau whose matrix is strictly lower triangular.
"""
struct ExplicitTableau <: TableauType end

"""
    abstract type ImplicitTableau <: TableauType end

Tableau whose matrix has at least one nonzero coefficient outside its strict
lower triangular part.
"""
abstract type ImplicitTableau <: TableauType end

"""
    struct DiagonallyImplicitTableau <: ImplicitTableau end

Tableau whose matrix is lower triangular, with at least one nonzero diagonal
coefficient.
"""
struct DiagonallyImplicitTableau <: ImplicitTableau end

"""
    struct FullyImplicitTableau <: ImplicitTableau end

Tableau whose matrix has at least one nonzero coefficient in its strict upper
triangular part.
"""
struct FullyImplicitTableau <: ImplicitTableau end

"""
    struct ImplicitExplicitTableau <: ImplicitTableau end

Pair of implicit and explicit tableaus that form a valid implicit-explicit
scheme.
"""
struct ImplicitExplicitTableau <: TableauType end

###################
# AbstractTableau #
###################
"""
    abstract type AbstractTableau{T} <: GridapType end

Type that stores the Butcher tableau corresponding to a Runge-Kutta scheme.
"""
abstract type AbstractTableau{T<:TableauType} <: GridapType end

"""
    TableauType(::AbstractTableau) -> TableauType

Return the `TableauType` of the tableau.
"""
TableauType(::AbstractTableau{T}) where {T} = T

"""
    get_matrix(tableau::AbstractTableau) -> AbstractMatrix

Return the matrix of the tableau.
"""
function Algebra.get_matrix(tableau::AbstractTableau)
  @abstractmethod
end

"""
    get_weights(tableau::AbstractTableau) -> AbstractVector

Return the weights of the tableau.
"""
function ReferenceFEs.get_weights(tableau::AbstractTableau)
  @abstractmethod
end

"""
    get_nodes(tableau::AbstractTableau) -> AbstractVector

Return the nodes of the tableau.
"""
function ReferenceFEs.get_nodes(tableau::AbstractTableau)
  @abstractmethod
end

"""
    get_order(tableau::AbstractTableau) -> Integer

Return the order of the scheme corresponding to the tableau.
"""
function Polynomials.get_order(tableau::AbstractTableau)
  @abstractmethod
end

##################
# GenericTableau #
##################
"""
    struct GenericTableau <: AbstractTableau end

Generic type that stores any type of Butcher tableau.
"""
struct GenericTableau{T<:TableauType} <: AbstractTableau{T}
  matrix::Matrix
  weights::Vector
  nodes::Vector
  order::Integer

  function GenericTableau(matrix, weights, order)
    nodes = reshape(sum(matrix, dims=2), size(matrix, 1))
    T = _tableau_type(matrix)
    new{T}(matrix, weights, nodes, order)
  end
end

function Algebra.get_matrix(tableau::GenericTableau)
  tableau.matrix
end

function ReferenceFEs.get_weights(tableau::GenericTableau)
  tableau.weights
end

function ReferenceFEs.get_nodes(tableau::GenericTableau)
  tableau.nodes
end

function Polynomials.get_order(tableau::GenericTableau)
  tableau.order
end

function _tableau_type(matrix::Matrix)
  T = ExplicitTableau
  n = size(matrix, 1)
  for i in 1:n
    if (i < n) && !iszero(matrix[i, i+1])
      T = FullyImplicitTableau
      break
    elseif !iszero(matrix[i, i])
      T = DiagonallyImplicitTableau
    end
  end
  T
end

###################
# EmbeddedTableau #
###################
"""
    struct EmbeddedTableau <: AbstractTableau end

Generic type that stores any type of embedded Butcher tableau.
"""
struct EmbeddedTableau{T} <: AbstractTableau{T}
  tableau::AbstractTableau{T}
  emb_weights::Vector
  emb_order::Integer
end

function Algebra.get_matrix(tableau::EmbeddedTableau)
  get_matrix(tableau.tableau)
end

function ReferenceFEs.get_weights(tableau::EmbeddedTableau)
  get_weights(tableau.tableau)
end

function ReferenceFEs.get_nodes(tableau::EmbeddedTableau)
  get_nodes(tableau.tableau)
end

function Polynomials.get_order(tableau::EmbeddedTableau)
  get_order(tableau.tableau)
end

"""
    get_embedded_weights(tableau::EmbeddedTableau) -> AbstractVector

Return the embedded weight of the tableau.
"""
function get_embedded_weights(tableau::EmbeddedTableau)
  tableau.emb_weights
end

"""
    get_embedded_order(tableau::EmbeddedTableau) -> Integer

Return the embedded order of the tableau.
"""
function get_embedded_order(tableau::EmbeddedTableau)
  tableau.emb_order
end

###############
# IMEXTableau #
###############
"""
    struct IMEXTableau <: AbstractTableau end

Generic type that stores any type of implicit-explicit pair of Butcher tableaus,
that form a valid IMEX scheme.
"""
struct IMEXTableau <: AbstractTableau{ImplicitExplicitTableau}
  im_tableau::AbstractTableau{<:ImplicitTableau}
  ex_tableau::AbstractTableau{ExplicitTableau}
  imex_order::Integer

  function IMEXTableau(im_tableau, ex_tableau, imex_order)
    Tim = TableauType(im_tableau)
    Tex = TableauType(ex_tableau)

    msg = """Invalid IMEX tableau:
    the first tableau must be implicit and the second must be explicit."""
    @assert (Tim <: ImplicitTableau && Tex == ExplicitTableau) msg

    msg = """Invalid IMEX tableau:
    the nodes of the implicit and explicit tableaus must coincide."""
    @assert isapprox(get_nodes(im_tableau), get_nodes(ex_tableau)) msg

    new(im_tableau, ex_tableau, imex_order)
  end
end

function Polynomials.get_order(tableau::IMEXTableau)
  tableau.imex_order
end

function get_imex_tableaus(tableau::IMEXTableau)
  (tableau.im_tableau, tableau.ex_tableau)
end

############################
# Concrete implementations #
############################
"""
    abstract type TableauName <: GridapType end

Name of a Butcher tableau.
"""
abstract type TableauName <: GridapType end

"""
    ButcherTableau(name::TableauName, type::Type) -> AbtractTableau

Builds the Butcher tableau corresponding to a `TableauName`.
"""
function ButcherTableau(name::TableauName, type::Type=Float64)
  @abstractmethod
end

function ButcherTableau(name::Symbol, type::Type=Float64)
  eval(:(ButcherTableau($name(), $type)))
end

# Redirect to ButcherTableau
for tableautype in (:GenericTableau, :EmbeddedTableau, :IMEXTableau)
  @eval begin
    function ($tableautype)(name::TableauName, type::Type=Float64)
      ButcherTableau(name, type)
    end
    function ($tableautype)(name::Symbol, type::Type=Float64)
      eval(:(ButcherTableau($name(), $type)))
    end
  end
end

##################################################
# Families of DIRK schemes with order conditions #
##################################################
function DIRK11(α::Real, ::Type{T}=Float64) where {T<:Real}
  matrix = T[α;;]
  weights = T[1]
  cond2 = (α ≈ 1 / 2)
  order = cond2 ? 2 : 1
  GenericTableau(matrix, weights, order)
end

function DIRK12(T::Type{<:Real}=Float64)
  RK11(1 / 2, T)
end

function DIRK22(α::Real, β::Real, γ::Real, ::Type{T}=Float64) where {T<:Real}
  δ = β - γ
  θ = (1 - 2 * α) / 2 / (β - α)
  matrix = T[
    α 0
    δ γ
  ]
  weights = T[1-θ, θ]
  cond31 = ((1 - θ) * α^2 + θ * β^2 ≈ 1 / 3)
  cond32 = ((1 - θ) * α^2 + θ * (δ * α + γ * β) ≈ 1 / 6)
  cond3 = cond31 && cond32
  order = cond3 ? 3 : 2
  GenericTableau(matrix, weights, order)
end

function DIRK23(λ::Real, ::Type{T}=Float64) where {T<:Real}
  α = 1 / 2 - sqrt(3) / 6 / λ
  β = sqrt(3) / 3 * λ
  γ = 1 / 2 - sqrt(3) / 6 * λ
  θ = 1 / (λ^2 + 1)
  matrix = T[
    α 0
    β γ
  ]
  weights = T[1-θ, θ]
  order = 3
  GenericTableau(matrix, weights, order)
end

####################
# Explicit schemes #
####################
"""
Forward Euler

Type              Explicit
Number of stages  1
Order             1
Stage order       1
"""
struct FE_1_0_1 <: TableauName end

function ButcherTableau(::FE_1_0_1, T::Type{<:Real}=Float64)
  DIRK11(0, T)
end

"""
3rd-order Strong Stability Preserving Runge-Kutta
SSPRK33

Type              Explicit
Number of stages  3
Order             3
Stage order       1
"""
struct SSPRK_3_0_3 <: TableauName end

function ButcherTableau(::SSPRK_3_0_3, ::Type{T}=Float64) where {T<:Real}
  a = 1
  b = 1 / 4
  c = 1 / 6
  d = 2 / 3
  matrix = T[
    0 0 0
    a 0 0
    b b 0
  ]
  weights = T[c, c, d]
  order = 3
  GenericTableau(matrix, weights, order)
end

###############################
# Diagonally-Implicit schemes #
###############################
"""
Backward Euler

Type              Diagonally Implicit
Number of stages  1
Order             1
Stage order       1
"""
struct BE_1_0_1 <: TableauName end

function ButcherTableau(::BE_1_0_1, T::Type{<:Real}=Float64)
  DIRK11(1, T)
end

"""
Crank-Nicolson
Trapezoidal rule
Lobatto IIIA2

Type              Diagonally Implicit
Number of stages  2
Order             2
Stage order       2
"""
struct CN_2_0_2 <: TableauName end

function ButcherTableau(::CN_2_0_2, T::Type{<:Real}=Float64)
  DIRK22(0, 1 / 2, 1 / 2, T)
end

"""
Qin and Zhang's SDIRK

Type              Singly Diagonally Implicit (Symplectic)
Number of stages  2
Order             2
Stage order       2
"""
struct SDIRK_2_0_2 <: TableauName end

function ButcherTableau(::SDIRK_2_0_2, T::Type{<:Real}=Float64)
  DIRK22(1 / 4, 3 / 4, 1 / 4, T)
end

"""
3rd order SDIRK

Type              Singly Diagonally Implicit
Number of stages  2
Order             3
Stage order       2
"""
struct SDIRK_2_0_3 <: TableauName end

function ButcherTableau(::SDIRK_2_0_3, T::Type{<:Real}=Float64)
  DIRK23(1, T)
end

"""
3rd order ESDIRK

Type              Diagonally Implicit
Number of stages  3
Order             2
Stage order       1
Embedded order    1
"""
struct ESDIRK_3_1_2 <: TableauName end

function ButcherTableau(::ESDIRK_3_1_2, ::Type{T}=Float64) where {T<:Real}
  c = (2 - sqrt(2)) / 2
  b = (1 - 2 * c) / (4 * c)
  a = 1 - b - c
  matrix = T[
    0 0 0
    c c 0
    a b c
  ]
  weights = T[a, b, c]
  order = 2
  tableau = GenericTableau(matrix, weights, order)

  ĉ = -2 * (c^2) * (1 - c + c^2) / (2 * c - 1)
  b̂ = c * (-2 + 7 * c - 5(c^2) + 4(c^3)) / (2 * (2 * c - 1))
  â = 1 - b̂ - ĉ
  emb_weights = T[â, b̂, ĉ]
  emb_order = 1
  EmbeddedTableau(tableau, emb_weights, emb_order)
end

"""
Trapezoidal Rule with Second Order Backward Difference Formula
TR-BDF2

Type              Diagonally Implicit
Number of stages  3
Order             3
Stage order       2
Embedded order    2
"""
struct TRBDF2_3_2_3 <: TableauName end

function ButcherTableau(::TRBDF2_3_2_3, ::Type{T}=Float64) where {T<:Real}
  γ = 2 - sqrt(2)
  d = γ / 2
  w = sqrt(2) / 4
  matrix = T[
    0 0 0
    d d 0
    w w d
  ]
  weights = T[w, w, d]
  order = 3
  tableau = GenericTableau(matrix, weights, order)

  ĉ = d / 3
  b̂ = (3 * w + 1) / 3
  â = 1 - b̂ - ĉ
  emb_weights = [â, b̂, ĉ]
  emb_order = 2
  EmbeddedTableau(tableau, emb_weights, emb_order)
end

"""
Double Trapezoidal Rule
TR-X2

Type              Diagonally Implicit
Number of stages  3
Order             3
Stage order       2
Embedded order    2
"""
struct TRX2_3_2_3 <: TableauName end

function ButcherTableau(::TRX2_3_2_3, ::Type{T}=Float64) where {T<:Real}
  a = 1 / 4
  b = 1 / 2
  matrix = T[
    0 0 0
    a a 0
    a b a
  ]
  weights = T[a, b, a]
  order = 3
  tableau = GenericTableau(matrix, weights, order)

  c = 1 / 6
  d = 2 / 3
  emb_weights = [c, d, c]
  emb_order = 2
  EmbeddedTableau(tableau, emb_weights, emb_order)
end

#############################
# Implicit-Explicit schemes #
#############################
"""
IMEX Forward-Backward-Euler

Type              Implicit-Explicit
Number of stages  2
Order             1
Stage order       1
"""
struct IMEX_FE_BE_2_0_1 <: TableauName end

function ButcherTableau(::IMEX_FE_BE_2_0_1, T::Type{<:Real}=Float64)
  im_tableau = DIRK22(0, 1, 1, T)
  ex_tableau = DIRK22(0, 1, 0, T)
  imex_order = 1
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end

"""
IMEX Midpoint

Type              Implicit-Explicit
Number of stages  2
Order             2
Stage order       2
"""
struct IMEX_Midpoint_2_0_2 <: TableauName end

function ButcherTableau(::IMEX_Midpoint_2_0_2, T::Type{<:Real}=Float64)
  im_tableau = DIRK22(0, 1 / 2, 1 / 2, T)
  ex_tableau = DIRK22(0, 1 / 2, 0, T)
  imex_order = 2
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end

const available_tableaus = [
  :FE_1_0_1
  :SSPRK_3_0_3
  :BE_1_0_1
  :CN_2_0_2
  :SDIRK_2_0_2
  :SDIRK_2_0_3
  :ESDIRK_3_1_2
  :TRBDF2_3_2_3
  :TRX2_3_2_3
  :IMEX_FE_BE_2_0_1
  :IMEX_Midpoint_2_0_2
]
