##########
# RKType #
##########
abstract type RKType <: GridapType end

struct ERK <: RKType end

abstract type IRK <: RKType end
struct FIRK <: IRK end
struct DIRK <: IRK end

struct IMEXRK <: RKType end

##########################
# AbstractButcherTableau #
##########################
abstract type AbstractButcherTableau{T<:RKType} <: GridapType end

RKType(::AbstractButcherTableau{T}) where {T} = T

function get_matrix(tableau::AbstractButcherTableau)
  @abstractmethod
end

function get_weights(tableau::AbstractButcherTableau)
  @abstractmethod
end

function get_nodes(tableau::AbstractButcherTableau)
  @abstractmethod
end

function get_order(tableau::AbstractButcherTableau)
  @abstractmethod
end

##################
# ButcherTableau #
##################
struct ButcherTableau{T} <: AbstractButcherTableau{T}
  matrix::Matrix
  weights::Vector
  nodes::Vector
  order::Int

  function ButcherTableau(matrix, weights, order)
    nodes = reshape(sum(matrix, dims=2), size(matrix, 1))
    T = _butcher_tableau_type(matrix)
    new{T}(matrix, weights, nodes, order)
  end
end

function get_matrix(tableau::ButcherTableau)
  tableau.matrix
end

function get_weights(tableau::ButcherTableau)
  tableau.weights
end

function get_nodes(tableau::ButcherTableau)
  tableau.nodes
end

function get_order(tableau::ButcherTableau)
  tableau.order
end

function _butcher_tableau_type(matrix::Matrix)
  T = ERK
  n = size(matrix, 1)
  for i in 1:n
    if (i < n) && !iszero(matrix[i, i+1])
      T = FIRK
      break
    elseif !iszero(matrix[i, i])
      T = DIRK
    end
  end
  T
end

##########################
# EmbeddedButcherTableau #
##########################
struct EmbeddedButcherTableau{T} <: AbstractButcherTableau{T}
  tableau::AbstractButcherTableau{T}
  emb_weights::Vector
  emb_order::Int
end

function get_matrix(tableau::EmbeddedButcherTableau)
  get_matrix(tableau.tableau)
end

function get_weights(tableau::EmbeddedButcherTableau)
  get_weights(tableau.tableau)
end

function get_nodes(tableau::EmbeddedButcherTableau)
  get_nodes(tableau.tableau)
end

function get_order(tableau::EmbeddedButcherTableau)
  get_order(tableau.tableau)
end

function get_embedded_weights(tableau::EmbeddedButcherTableau)
  tableau.emb_weights
end

function get_embedded_order(tableau::EmbeddedButcherTableau)
  tableau.emb_order
end

######################
# IMEXButcherTableau #
######################
struct IMEXButcherTableau <: AbstractButcherTableau{IMEXRK}
  im_tableau::AbstractButcherTableau
  ex_tableau::AbstractButcherTableau
  order::Int

  function IMEXButcherTableau(im_tableau, ex_tableau, order)
    Tim = RKType(im_tableau)
    Tex = RKType(ex_tableau)

    msg = """Invalid IMEX Butcher tableau:
    the first tableau must be implicit and the second must be explicit."""
    @assert (Tim <: IRK && Tex == ERK) msg

    msg = """Invalid IMEX Butcher tableau:
    the nodes of the implicit and explicit tableaus must coincide."""
    @assert isapprox(get_nodes(im_tableau), get_nodes(ex_tableau)) msg

    new(im_tableau, ex_tableau, order)
  end
end

function get_matrix(tableau::IMEXButcherTableau)
  get_matrix(tableau.im_tableau), get_matrix(tableau.ex_tableau)
end

function get_weights(tableau::IMEXButcherTableau)
  get_weights(tableau.im_tableau), get_weights(tableau.ex_tableau)
end

function get_nodes(tableau::IMEXButcherTableau)
  get_nodes(tableau.im_tableau)
end

function get_order(tableau::IMEXButcherTableau)
  tableau.order
end

############################
# Concrete implementations #
############################
abstract type ButcherTableauName <: GridapType end

# Redirect to ButcherTableau
function ButcherTableau(name::Symbol, type::Type=Float64)
  eval(:(ButcherTableau($name(), $type)))
end

function EmbeddedButcherTableau(name::Symbol, type::Type=Float64)
  eval(:(ButcherTableau($name(), $type)))
end

function IMEXButcherTableau(name::Symbol, type::Type=Float64)
  eval(:(ButcherTableau($name(), $type)))
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
struct FE_1_0_1 <: ButcherTableauName end

function ButcherTableau(::FE_1_0_1, ::Type{T}=Float64) where {T}
  matrix = T[0;;]
  weights = T[1]
  order = 1
  ButcherTableau(matrix, weights, order)
end

"""
3rd-order Strong Stability Preserving Runge-Kutta
SSPRK33

Type              Explicit
Number of stages  3
Order             3
Stage order       1
"""
struct SSPRK_3_0_3 <: ButcherTableauName end

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
  ButcherTableau(matrix, weights, order)
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
struct BE_1_0_1 <: ButcherTableauName end

function ButcherTableau(::BE_1_0_1, ::Type{T}=Float64) where {T}
  matrix = T[1;;]
  weights = T[1]
  order = 1
  ButcherTableau(matrix, weights, order)
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
struct CN_2_0_2 <: ButcherTableauName end

function ButcherTableau(::CN_2_0_2, ::Type{T}=Float64) where {T}
  a = 1 / 2
  matrix = T[
    0 0
    a a
  ]
  weights = T[a, a]
  order = 2
  ButcherTableau(matrix, weights, order)
end

"""
Qin and Zhang's SDIRK

Type              Singly Diagonally Implicit (Symplectic)
Number of stages  2
Order             2
Stage order       2
"""
struct SDIRK_2_0_2 <: ButcherTableauName end

function ButcherTableau(::SDIRK_2_0_2, ::Type{T}=Float64) where {T}
  a = 1 / 4
  b = 1 / 2
  matrix = T[
    a 0
    b a
  ]
  weights = T[b, b]
  order = 2
  ButcherTableau(matrix, weights, order)
end

"""
3rd order SDIRK

Type              Diagonally Implicit
Number of stages  2
Order             3
Stage order       2
"""
struct SDIRK_2_0_3 <: ButcherTableauName end

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
  ButcherTableau(matrix, weights, order)
end

"""
3rd order ESDIRK

Type              Diagonally Implicit
Number of stages  3
Order             2
Stage order       1
Embedded order    1
"""
struct ESDIRK_3_1_2 <: ButcherTableauName end

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
  tableau = ButcherTableau(matrix, weights, order)

  ĉ = -2 * (c^2) * (1 - c + c^2) / (2 * c - 1)
  b̂ = c * (-2 + 7 * c - 5(c^2) + 4(c^3)) / (2 * (2 * c - 1))
  â = 1 - b̂ - ĉ
  emb_weights = T[â, b̂, ĉ]
  emb_order = 1
  EmbeddedButcherTableau(tableau, emb_weights, emb_order)
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
struct TRBDF2_3_2_3 <: ButcherTableauName end

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
  tableau = ButcherTableau(matrix, weights, order)

  ĉ = d / 3
  b̂ = (3 * w + 1) / 3
  â = 1 - b̂ - ĉ
  emb_weights = [â, b̂, ĉ]
  emb_order = 2
  EmbeddedButcherTableau(tableau, emb_weights, emb_order)
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
struct TRX2_3_2_3 <: ButcherTableauName end

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
  tableau = ButcherTableau(matrix, weights, order)

  c = 1 / 6
  d = 2 / 3
  emb_weights = [c, d, c]
  emb_order = 2
  EmbeddedButcherTableau(tableau, emb_weights, emb_order)
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
struct IMEX_FE_BE_2_0_1 <: ButcherTableauName end

function ButcherTableau(::IMEX_FE_BE_2_0_1, ::Type{T}=Float64) where {T}
  a = 1
  im_matrix = T[
    0 0
    0 a
  ]
  im_weights = T[0, a]
  im_order = 1
  im_tableau = ButcherTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0
    a 0
  ]
  ex_weights = T[0, a]
  ex_order = 1
  ex_tableau = ButcherTableau(ex_matrix, ex_weights, ex_order)

  order = 1
  IMEXButcherTableau(im_tableau, ex_tableau, order)
end

"""
IMEX Midpoint

Type              Implicit-Explicit
Number of stages  2
Order             2
Stage order       2
"""
struct IMEX_Midpoint_2_0_2 <: ButcherTableauName end

function ButcherTableau(::IMEX_Midpoint_2_0_2, ::Type{T}=Float64) where {T}
  a = 1
  b = 1 / 2
  im_matrix = T[
    0 0
    0 b
  ]
  im_weights = T[0, a]
  im_order = 1
  im_tableau = ButcherTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0
    b 0
  ]
  ex_weights = T[0, a]
  ex_order = 1
  ex_tableau = ButcherTableau(ex_matrix, ex_weights, ex_order)

  order = 2
  IMEXButcherTableau(im_tableau, ex_tableau, order)
end
