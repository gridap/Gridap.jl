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
    if any(j -> !iszero(matrix[i, j]), i+1:n)
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
  is_padded::Bool

  function IMEXTableau(im_tableau, ex_tableau, imex_order)
    Tim = TableauType(im_tableau)
    Tex = TableauType(ex_tableau)

    msg = """Invalid IMEX tableau:
    the first tableau must be implicit and the second must be explicit."""
    @check (Tim <: ImplicitTableau && Tex == ExplicitTableau) msg

    msg = """Invalid IMEX tableau:
    the nodes of the implicit and explicit tableaus must coincide."""
    @check isapprox(get_nodes(im_tableau), get_nodes(ex_tableau)) msg

    is_padded = _is_padded(im_tableau)

    new(im_tableau, ex_tableau, imex_order, is_padded)
  end
end

function Polynomials.get_order(tableau::IMEXTableau)
  tableau.imex_order
end

"""
    get_imex_tableaus(tableau::IMEXTableau)

Return the pair of implicit and explicit tableaus of the given IMEX tableau.
"""
function get_imex_tableaus(tableau::IMEXTableau)
  (tableau.im_tableau, tableau.ex_tableau)
end

"""
    is_padded(tableau::IMEXTableau) -> Bool
"""
function is_padded(tableau::IMEXTableau)
  tableau.is_padded
end

function _is_padded(tableau::AbstractTableau)
  A = get_matrix(tableau)
  b = get_matrix(tableau)
  iszero(b[1]) && all(i -> iszero(A[i, 1]), axes(A, 1))
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

##################
# Import schemes #
##################

include("TableausEX.jl")

include("TableausDIM.jl")

include("TableausIMEX.jl")

"""
List of available Butcher tableaus.
"""
const available_tableaus = [
  :EXRK_Euler_1_1,
  :EXRK_Midpoint_2_2,
  :EXRK_SSP_2_2,
  :EXRK_Heun_2_2,
  :EXRK_Ralston_2_2,
  :EXRK_Kutta_3_3,
  :EXRK_Heun_3_3,
  :EXRK_Wray_3_3,
  :EXRK_VanDerHouwen_3_3,
  :EXRK_Ralston_3_3,
  :EXRK_SSP_3_3,
  :EXRK_SSP_3_2,
  :EXRK_Fehlberg_3_2,
  :EXRK_RungeKutta_4_4,
  :EXRK_Simpson_4_4,
  :EXRK_Ralston_4_4,
  :EXRK_SSP_4_3,
  :EXRK_BogackiShampine_4_3,
  :SDIRK_Euler_1_1,
  :SDIRK_Midpoint_1_2,
  :DIRK_CrankNicolson_2_2,
  :SDIRK_QinZhang_2_2,
  :DIRK_LobattoIIIA_2_2,
  :DIRK_RadauI_2_3,
  :DIRK_RadauII_2_3,
  :SDIRK_LobattoIIIC_2_2,
  :SDIRK_2_2,
  :SDIRK_SSP_2_3,
  :SDIRK_Crouzeix_2_3,
  :SDIRK_3_2,
  :DIRK_TRBDF_3_2,
  :DIRK_TRX_3_2,
  :SDIRK_3_3,
  :SDIRK_Crouzeix_3_4,
  :SDIRK_Norsett_3_4,
  :DIRK_LobattoIIIC_3_4,
  :SDIRK_4_3,
]

"""
List of available tableaus for IMEX schemes.
"""
const available_imex_tableaus = [
  :IMEXRK_1_1_1,
  :IMEXRK_1_2_1,
  :IMEXRK_1_2_2,
  :IMEXRK_2_2_2,
  :IMEXRK_2_3_2,
  :IMEXRK_2_3_3,
  :IMEXRK_3_4_3,
  :IMEXRK_4_4_3,
]
