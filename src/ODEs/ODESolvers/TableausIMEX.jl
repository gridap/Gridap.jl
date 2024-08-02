# All these schemes come from the following paper
# Implicit-explicit Runge-Kutta methods for time-dependent partial differential
# equations, Uri M. Ascher, Steven J. Ruuth, Raymond J. Spiteri, Applied
# numerical mathematics 1997.
# https://www.sciencedirect.com/science/article/abs/pii/S0168927497000561

############
# 2 stages #
############
"""
    IMEXRK_1_1_1

Backward-Forward Euler pair, order 1
"""
struct IMEXRK_1_1_1 <: TableauName end

function ButcherTableau(::IMEXRK_1_1_1, T::Type{<:Real}=Float64)
  im_matrix = T[
    0 0
    0 1
  ]
  im_weights = T[0, 1]
  im_order = 1
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0
    1 0
  ]
  ex_weights = T[1, 0]
  ex_order = 1
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  imex_order = 1
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end

"""
    IMEXRK_1_2_1

Backward-Forward Euler pair with same weights, order 1
"""
struct IMEXRK_1_2_1 <: TableauName end

function ButcherTableau(::IMEXRK_1_2_1, T::Type{<:Real}=Float64)
  im_matrix = T[
    0 0
    0 1
  ]
  im_weights = T[0, 1]
  im_order = 1
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0
    1 0
  ]
  ex_weights = T[0, 1]
  ex_order = 1
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  imex_order = 1
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end

"""
    IMEXRK_1_2_2

Implicit-Explicit midpoint pair, order 2
"""
struct IMEXRK_1_2_2 <: TableauName end

function ButcherTableau(::IMEXRK_1_2_2, T::Type{<:Real}=Float64)
  a = 1 / 2
  im_matrix = T[
    0 0
    0 a
  ]
  im_weights = T[0, 1]
  im_order = 2
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0
    a 0
  ]
  ex_weights = T[0, 1]
  ex_order = 2
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  imex_order = 2
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end

"""
    IMEXRK_2_2_2

L-stable, 2-stage, 2-order SDIRK
"""
struct IMEXRK_2_2_2 <: TableauName end

function ButcherTableau(::IMEXRK_2_2_2, T::Type{<:Real}=Float64)
  a = (2 - sqrt(2)) / 2
  b = 1 - a
  c = 1 - 1 / 2 / a
  d = 1 - c
  im_matrix = T[
    0 0 0
    0 a 0
    0 b a
  ]
  im_weights = T[0, b, a]
  im_order = 2
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0 0
    a 0 0
    c d 0
  ]
  ex_weights = T[c, d, 0]
  ex_order = 2
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  imex_order = 2
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end

"""
    IMEXRK_2_3_2

L-stable, 2-stage, 2-order SDIRK
"""
struct IMEXRK_2_3_2 <: TableauName end

function ButcherTableau(::IMEXRK_2_3_2, T::Type{<:Real}=Float64)
  a = (2 - sqrt(2)) / 2
  b = 1 - a
  c = -2 * sqrt(2) / 3
  d = 1 - c
  im_matrix = T[
    0 0 0
    0 a 0
    0 b a
  ]
  im_weights = T[0, b, a]
  im_order = 3
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0 0
    a 0 0
    c d 0
  ]
  ex_weights = T[0, b, a]
  ex_order = 3
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  imex_order = 3
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end

"""
    IMEXRK_2_3_3

2-stage, 3-order SDIRK scheme with best damping properties
"""
struct IMEXRK_2_3_3 <: TableauName end

function ButcherTableau(::IMEXRK_2_3_3, T::Type{<:Real}=Float64)
  a = (3 + sqrt(3)) / 6
  b = 1 - 2 * a
  c = 1 // 2
  d = a - 1
  e = 2 * (1 - a)
  im_matrix = T[
    0 0 0
    0 a 0
    0 b a
  ]
  im_weights = T[0, c, c]
  im_order = 3
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0 0
    a 0 0
    d e 0
  ]
  ex_weights = T[0, c, c]
  ex_order = 3
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  imex_order = 3
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end

"""
    IMEXRK_3_4_3

L-stable, 3-stage, 3-order SDIRK
"""
struct IMEXRK_3_4_3 <: TableauName end

function ButcherTableau(::IMEXRK_3_4_3, T::Type{<:Real}=Float64)
  # a is the middle root of 6 * x^3 - 18 * x^2 + 9 * x - 1
  a = 0.435866521508459
  b = (1 - a) / 2
  c = -3 * a^2 / 2 + 4 * a - 1 // 4
  d = 3 * a^2 / 2 - 5 * a + 5 / 4
  h = 0.5529291479
  e = (1 - 9 * a / 2 + 3 * a^2 / 2 + 11 // 4 - 21 * a / 2 + 15 * a^2 / 4) * h - 7 // 2 + 13 * a - 9 * a^2 / 2
  f = (-1 + 9 * a / 2 - 3 * a^2 / 2 - 11 // 4 + 21 * a / 2 - 15 * a^2 / 4) * h + 4 - 25 * a / 2 + 9 * a^2 / 2
  g = 1 - 2 * h
  im_matrix = T[
    0 0 0 0
    0 a 0 0
    0 b a 0
    0 c d a
  ]
  im_weights = T[0, c, d, a]
  im_order = 3
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0 0 0
    a 0 0 0
    e f 0 0
    g h h 0
  ]
  ex_weights = T[0, c, d, a]
  ex_order = 3
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  imex_order = 3
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end

"""
    IMEXRK_4_4_3

L-stable, 4-stage, 3-order SDIRK
"""
struct IMEXRK_4_4_3 <: TableauName end

function ButcherTableau(::IMEXRK_4_4_3, T::Type{<:Real}=Float64)
  a = 1 // 2
  b = 1 // 6
  c = -1 // 2
  d = 3 // 2
  e = -3 // 2
  f = 11 // 18
  g = 1 // 18
  h = 5 // 6
  i = -5 // 6
  j = 1 // 4
  k = 7 // 4
  l = 3 // 4
  m = -7 // 4
  im_matrix = T[
    0 0 0 0 0
    0 a 0 0 0
    0 b a 0 0
    0 c a a 0
    0 d e a a
  ]
  im_weights = T[0, d, e, a, a]
  im_order = 3
  im_tableau = GenericTableau(im_matrix, im_weights, im_order)

  ex_matrix = T[
    0 0 0 0 0
    a 0 0 0 0
    f g 0 0 0
    h i a 0 0
    j k l m 0
  ]
  ex_weights = T[j, k, l, m, 0]
  ex_order = 3
  ex_tableau = GenericTableau(ex_matrix, ex_weights, ex_order)

  imex_order = 3
  IMEXTableau(im_tableau, ex_tableau, imex_order)
end
