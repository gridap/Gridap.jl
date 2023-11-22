###############
# TableauType #
###############
"""
    abstract type TableauType <: GridapType end

Trait for `AbstractTableau` to differentiate between explicit, implicit and
implicit-explicit tableaus
"""
abstract type TableauType <: GridapType end

"""
    struct ExplicitTableau <: TableauType end

Tableau whose matrix is strictly lower triangular
"""
struct ExplicitTableau <: TableauType end

"""
    abstract type ImplicitTableau <: TableauType end

Tableau whose matrix has at least one nonzero coefficient outside the strict
lower triangular part
"""
abstract type ImplicitTableau <: TableauType end

"""
    struct DiagonallyImplicitTableau <: ImplicitTableau end

Tableau whose matrix is lower triangular, with at least one nonzero diagonal
coefficient
"""
struct DiagonallyImplicitTableau <: ImplicitTableau end

"""
    struct FullyImplicitTableau <: ImplicitTableau end

Tableau whose matrix has at least one nonzero coefficient in the strict upper
triangular part
"""
struct FullyImplicitTableau <: ImplicitTableau end

"""
    struct ImplicitExplicitTableau <: ImplicitTableau end

Pair of implicit and explicit tableaus that form a valid implicit-explicit
scheme
"""
struct ImplicitExplicitTableau <: TableauType end

###################
# AbstractTableau #
###################
"""
    abstract type AbstractTableau{T} <: GridapType end

Type that stores the Butcher tableau corresponding to a Runge-Kutta scheme
"""
abstract type AbstractTableau{T<:TableauType} <: GridapType end

"""
    TableauType(::AbstractTableau) -> TableauType

Return the `TableauType` of the `AbtractTableau`
"""
TableauType(::AbstractTableau{T}) where {T} = T

"""
    get_matrix(tableau::AbstractTableau) -> AbstractMatrix

Return the matrix of the tableau
"""
function Algebra.get_matrix(tableau::AbstractTableau)
  @abstractmethod
end

"""
    get_weights(tableau::AbstractTableau) -> AbstractVector

Return the weights of the tableau
"""
function ReferenceFEs.get_weights(tableau::AbstractTableau)
  @abstractmethod
end

"""
    get_nodes(tableau::AbstractTableau) -> AbstractVector

Return the nodes of the tableau
"""
function ReferenceFEs.get_nodes(tableau::AbstractTableau)
  @abstractmethod
end

"""
    get_order(tableau::AbstractTableau) -> Integer

Return the order of the scheme corresponding to the tableau
"""
function Polynomials.get_order(tableau::AbstractTableau)
  @abstractmethod
end

##################
# GenericTableau #
##################
"""
    struct GenericTableau <: AbstractTableau

Generic type that stores any type of Butcher tableau
"""
struct GenericTableau{T<:TableauType} <: AbstractTableau{T}
  matrix::Matrix
  weights::Vector
  nodes::Vector
  order::Int

  function GenericTableau(matrix, weights, order)
    nodes = reshape(sum(matrix, dims=2), size(matrix, 1))
    T = _butcher_tableau_type(matrix)
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

function _butcher_tableau_type(matrix::Matrix)
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
    struct EmbeddedTableau <: AbstractTableau

Generic type that stores any type of embedded Butcher tableau
"""
struct EmbeddedTableau{T} <: AbstractTableau{T}
  tableau::AbstractTableau{T}
  emb_weights::Vector
  emb_order::Int
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

Return the embedded weight of the tableau
"""
function get_embedded_weights(tableau::EmbeddedTableau)
  tableau.emb_weights
end

"""
    get_embedded_order(tableau::EmbeddedTableau) -> Integer

Return the embedded order of the tableau
"""
function get_embedded_order(tableau::EmbeddedTableau)
  tableau.emb_order
end

###############
# IMEXTableau #
###############
"""
    struct IMEXTableau <: AbstractTableau

Generic type that stores any type of implicit-explicit pair of Butcher tableaus,
that form a valid IMEX scheme
"""
struct IMEXTableau <: AbstractTableau{ImplicitExplicitTableau}
  im_tableau::AbstractTableau
  ex_tableau::AbstractTableau
  order::Int

  function IMEXTableau(im_tableau, ex_tableau, order)
    Tim = TableauType(im_tableau)
    Tex = TableauType(ex_tableau)

    msg = """Invalid IMEX tableau:
    the first tableau must be implicit and the second must be explicit."""
    @assert (Tim <: ImplicitTableau && Tex == ExplicitTableau) msg

    msg = """Invalid IMEX tableau:
    the nodes of the implicit and explicit tableaus must coincide."""
    @assert isapprox(get_nodes(im_tableau), get_nodes(ex_tableau)) msg

    new(im_tableau, ex_tableau, order)
  end
end

function Algebra.get_matrix(tableau::IMEXTableau)
  get_matrix(tableau.im_tableau), get_matrix(tableau.ex_tableau)
end

function ReferenceFEs.get_weights(tableau::IMEXTableau)
  get_weights(tableau.im_tableau), get_weights(tableau.ex_tableau)
end

function ReferenceFEs.get_nodes(tableau::IMEXTableau)
  get_nodes(tableau.im_tableau)
end

function Polynomials.get_order(tableau::IMEXTableau)
  tableau.order
end

############################
# Concrete implementations #
############################
"""
    abstract type TableauName <: GridapType end

Name of a Butcher tableau
"""
abstract type TableauName <: GridapType end

"""
    ButcherTableau(name::TableauName, type::Type) -> AbtractTableau

Builds the Butcher tableau corresponding to a `TableauName`
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

function ButcherTableau(::FE_1_0_1, ::Type{T}=Float64) where {T}
  matrix = T[0;;]
  weights = T[1]
  order = 1
  GenericTableau(matrix, weights, order)
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

function ButcherTableau(::SSPRK_3_0_3, ::Type{T}=Float64) where {T}
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

function ButcherTableau(::BE_1_0_1, ::Type{T}=Float64) where {T}
  matrix = T[1;;]
  weights = T[1]
  order = 1
  GenericTableau(matrix, weights, order)
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

function ButcherTableau(::CN_2_0_2, ::Type{T}=Float64) where {T}
  a = 1 / 2
  matrix = T[
    0 0
    a a
  ]
  weights = T[a, a]
  order = 2
  GenericTableau(matrix, weights, order)
end

"""
Qin and Zhang's SDIRK

Type              Singly Diagonally Implicit (Symplectic)
Number of stages  2
Order             2
Stage order       2
"""
struct SDIRK_2_0_2 <: TableauName end

function ButcherTableau(::SDIRK_2_0_2, ::Type{T}=Float64) where {T}
  a = 1 / 4
  b = 1 / 2
  matrix = T[
    a 0
    b a
  ]
  weights = T[b, b]
  order = 2
  GenericTableau(matrix, weights, order)
end

"""
3rd order SDIRK

Type              Diagonally Implicit
Number of stages  2
Order             3
Stage order       2
"""
struct SDIRK_2_0_3 <: TableauName end

function ButcherTableau(::SDIRK_2_0_3, ::Type{T}=Float64) where {T}
  a = (3 - sqrt(3)) / 6
  b = 1 - 2 * a
  c = 1 / 2
  matrix = T[
    a 0
    b a
  ]
  weights = T[c, c]
  order = 3
  GenericTableau(matrix, weights, order)
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

function ButcherTableau(::ESDIRK_3_1_2, ::Type{T}=Float64) where {T}
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

function ButcherTableau(::TRBDF2_3_2_3, ::Type{T}=Float64) where {T}
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

function ButcherTableau(::TRX2_3_2_3, ::Type{T}=Float64) where {T}
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

function ButcherTableau(::IMEX_FE_BE_2_0_1, ::Type{T}=Float64) where {T}
  a = 1
  im_matrix = T[
    0 0
    0 a
  ]
  im_weights = T[0, a]
  im_order = 1
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0
    a 0
  ]
  ex_weights = T[0, a]
  ex_order = 1
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  order = 1
  IMEXTableau(im_tableau, ex_tableau, order)
end

"""
IMEX Midpoint

Type              Implicit-Explicit
Number of stages  2
Order             2
Stage order       2
"""
struct IMEX_Midpoint_2_0_2 <: TableauName end

function ButcherTableau(::IMEX_Midpoint_2_0_2, ::Type{T}=Float64) where {T}
  a = 1
  b = 1 / 2
  im_matrix = T[
    0 0
    0 b
  ]
  im_weights = T[0, a]
  im_order = 1
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0
    b 0
  ]
  ex_weights = T[0, a]
  ex_order = 1
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  order = 2
  IMEXTableau(im_tableau, ex_tableau, order)
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
