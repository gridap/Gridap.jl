###########
# 1 stage #
###########
function EXRK11(::Type{T}) where {T}
  matrix = T[0;;]
  weights = T[1]
  order = 1
  GenericTableau(matrix, weights, order)
end

"""
    EXRK_Euler_1_1
    FE
"""
struct EXRK_Euler_1_1 <: TableauName end

function ButcherTableau(::EXRK_Euler_1_1, ::Type{T}=Float64) where {T}
  EXRK11(T)
end

############
# 2 stages #
############
function EXRK22(α::Real, ::Type{T}=Float64) where {T}
  a = 1 // 2 / α
  matrix = T[
    0 0
    α 0
  ]
  weights = T[1-a, a]
  order = 2
  GenericTableau(matrix, weights, order)
end

"""
    EXRK_Midpoint_2_2
"""
struct EXRK_Midpoint_2_2 <: TableauName end

function ButcherTableau(::EXRK_Midpoint_2_2, ::Type{T}=Float64) where {T}
  EXRK22(1 // 2, T)
end

"""
    EXRK_SSP_2_2
    EXRK_Heun_2_2
"""
struct EXRK_SSP_2_2 <: TableauName end

function ButcherTableau(::EXRK_SSP_2_2, ::Type{T}=Float64) where {T}
  tableau = EXRK22(1, T)
  emb_weights = T[1, 0]
  emb_order = 1
  EmbeddedTableau(tableau, emb_weights, emb_order)
end

struct EXRK_Heun_2_2 <: TableauName end

function ButcherTableau(::EXRK_Heun_2_2, ::Type{T}=Float64) where {T}
  ButcherTableau(EXRK_SSP_2_2(), T)
end

"""
    EXRK_Ralston_2_2
"""
struct EXRK_Ralston_2_2 <: TableauName end

function ButcherTableau(::EXRK_Ralston_2_2, ::Type{T}=Float64) where {T}
  EXRK22(2 // 3, T)
end

############
# 3 stages #
############
function EXRK33(α::Real, β::Real, ::Type{T}=Float64) where {T}
  b = β * (β - α) / α / (2 - 3 * α)
  c = β - b
  d = (3 * β - 2) // 6 / α / (β - α)
  e = (2 - 3 * α) // 6 / β / (β - α)
  matrix = T[
    0 0 0
    α 0 0
    c b 0
  ]
  weights = T[1-d-e, d, e]
  order = 3
  GenericTableau(matrix, weights, order)
end

function EXRK33_1(α::Real, ::Type{T}=Float64) where {T}
  a = 2 // 3
  b = 1 // 4 / α
  c = 2 // 3 - b
  d = 1 // 4
  matrix = T[
    0 0 0
    a 0 0
    c b 0
  ]
  weights = T[d, 1-α-d, α]
  order = 3
  GenericTableau(matrix, weights, order)
end

function EXRK33_2(α::Real, ::Type{T}=Float64) where {T}
  a = 2 // 3
  b = 1 // 4 / α
  c = -b
  d = 3 // 4
  matrix = T[
    0 0 0
    a 0 0
    c b 0
  ]
  weights = T[1-α-d, d, α]
  order = 3
  GenericTableau(matrix, weights, order)
end

"""
    EXRK_Kutta_3_3
"""
struct EXRK_Kutta_3_3 <: TableauName end

function ButcherTableau(::EXRK_Kutta_3_3, ::Type{T}=Float64) where {T}
  EXRK33(1 // 2, 1, T)
end

"""
    EXRK_Heun_3_3
"""
struct EXRK_Heun_3_3 <: TableauName end

function ButcherTableau(::EXRK_Heun_3_3, ::Type{T}=Float64) where {T}
  EXRK33(1 // 3, 2 // 3, T)
end

"""
    EXRK_Wray_3_3
    EXRK_VanDerHouwen_3_3
"""
struct EXRK_Wray_3_3 <: TableauName end

function ButcherTableau(::EXRK_Wray_3_3, ::Type{T}=Float64) where {T}
  EXRK33(8 // 15, 2 // 3, T)
end

struct EXRK_VanDerHouwen_3_3 <: TableauName end

function ButcherTableau(::EXRK_VanDerHouwen_3_3, ::Type{T}=Float64) where {T}
  ButcherTableau(EXRK_Wray_3_3(), T)
end

"""
    EXRK_Ralston_3_3
"""
struct EXRK_Ralston_3_3 <: TableauName end

function ButcherTableau(::EXRK_Ralston_3_3, ::Type{T}=Float64) where {T}
  EXRK33(1 // 2, 3 // 4, T)
end

"""
    EXRK_SSP_3_3
"""
struct EXRK_SSP_3_3 <: TableauName end

function ButcherTableau(::EXRK_SSP_3_3, ::Type{T}=Float64) where {T}
  EXRK33(1, 1 // 2, T)
end

"""
    EXRK_SSP_3_2
"""
struct EXRK_SSP_3_2 <: TableauName end

function ButcherTableau(::EXRK_SSP_3_2, ::Type{T}=Float64) where {T}
  a = 1 // 2
  b = 1 // 3
  matrix = T[
    0 0 0
    a 0 0
    a a 0
  ]
  weights = T[b, b, b]
  order = 2
  GenericTableau(matrix, weights, order)
end

"""
    EXRK_Fehlberg_3_2
"""
struct EXRK_Fehlberg_3_2 <: TableauName end

function ButcherTableau(::EXRK_Fehlberg_3_2, ::Type{T}=Float64) where {T}
  a = 1 // 2
  b = 1 // 256
  c = 255 // 256
  d = 1 // 512
  matrix = T[
    0 0 0
    a 0 0
    b c 0
  ]
  weights = T[d, c, d]
  order = 2
  tableau = GenericTableau(matrix, weights, order)

  emb_weights = T[b, c, 0]
  emb_order = 1
  EmbeddedTableau(tableau, emb_weights, emb_order)
end

############
# 4 stages #
############
"""
    EXRK_RungeKutta_4_4
"""
struct EXRK_RungeKutta_4_4 <: TableauName end

function ButcherTableau(::EXRK_RungeKutta_4_4, ::Type{T}=Float64) where {T}
  a = 1 // 2
  b = 1 // 6
  c = 1 // 3
  matrix = T[
    0 0 0 0
    a 0 0 0
    0 a 0 0
    0 0 1 0
  ]
  weights = T[b, c, c, b]
  order = 4
  GenericTableau(matrix, weights, order)
end

"""
    EXRK_Simpson_4_4
"""
struct EXRK_Simpson_4_4 <: TableauName end

function ButcherTableau(::EXRK_Simpson_4_4, ::Type{T}=Float64) where {T}
  a = 1 // 3
  b = -1 // 3
  c = -1
  d = 1 // 8
  e = 3 // 8
  matrix = T[
    0 0 0 0
    a 0 0 0
    b 1 0 0
    1 c 1 0
  ]
  weights = T[d, e, e, d]
  order = 4
  GenericTableau(matrix, weights, order)
end

"""
    EXRK_Ralston_4_4
"""
struct EXRK_Ralston_4_4 <: TableauName end

function ButcherTableau(::EXRK_Ralston_4_4, ::Type{T}=Float64) where {T}
  α = 2 // 5
  β = 7 // 8 - 3 * sqrt(5) / 16
  a = β * (β - α) / 2 / α / (1 - 2 * α)
  b = β - a
  c = (1 - α) * (α + β - 1 - (2 * β - 1)^2) / 2 / α / (β - α) / (6 * α * β - 4 * (α + β) + 3)
  d = (1 - 2 * α) * (1 - α) * (1 - β) / β / (β - α) / (6 * α * β - 4 * (α + β) + 3)
  e = 1 - c - d
  f = 1 // 2 + (1 - 2 * (α + β)) / 12 / α / β
  g = (2 * β - 1) / 12 / α / (β - α) / (1 - α)
  h = (1 - 2 * α) / 12 / β / (β - α) / (1 - β)
  i = 1 // 2 + (2 * (α + β) - 3) / 12 / (1 - α) / (1 - β)
  matrix = T[
    0 0 0 0
    α 0 0 0
    b a 0 0
    e c d 0
  ]
  weights = T[f, g, h, i]
  order = 4
  GenericTableau(matrix, weights, order)
end

"""
    EXRK_SSP_4_3
"""
struct EXRK_SSP_4_3 <: TableauName end

function ButcherTableau(::EXRK_SSP_4_3, ::Type{T}=Float64) where {T}
  a = 1 // 2
  b = 1 // 6
  c = 1 // 4
  matrix = T[
    0 0 0 0
    a 0 0 0
    a a 0 0
    b b b 0
  ]
  weights = T[b, b, b, a]
  order = 3
  tableau = GenericTableau(matrix, weights, order)

  emb_weights = T[c, c, c, c]
  emb_order = 2
  EmbeddedTableau(tableau, emb_weights, emb_order)
end

"""
    EXRK_BogackiShampine_4_3
"""
struct EXRK_BogackiShampine_4_3 <: TableauName end

function ButcherTableau(::EXRK_BogackiShampine_4_3, ::Type{T}=Float64) where {T}
  a = 1 // 2
  b = 3 // 4
  c = 2 // 9
  d = 1 // 3
  e = 4 // 9
  matrix = T[
    0 0 0 0
    a 0 0 0
    0 b 0 0
    c d e 0
  ]
  weights = T[c, d, e, 0]
  order = 3
  tableau = GenericTableau(matrix, weights, order)

  emb_weights = T[0, a, b, 1]
  emb_order = 2
  EmbeddedTableau(tableau, emb_weights, emb_order)
end
