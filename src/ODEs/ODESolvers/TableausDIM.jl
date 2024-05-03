###########
# 1 stage #
###########
function SDIRK11(α::Real, ::Type{T}=Float64) where {T<:Real}
  matrix = T[α;;]
  weights = T[1]
  cond2 = (α ≈ 1 / 2)
  order = cond2 ? 2 : 1
  GenericTableau(matrix, weights, order)
end

function SDIRK12(T::Type{<:Real}=Float64)
  SDIRK11(1 / 2, T)
end

"""
    SDIRK_Euler_1_1
"""
struct SDIRK_Euler_1_1 <: TableauName end

function ButcherTableau(::SDIRK_Euler_1_1, ::Type{T}=Float64) where {T}
  SDIRK11(1, T)
end

"""
    SDIRK_Midpoint_1_2
"""
struct SDIRK_Midpoint_1_2 <: TableauName end

function ButcherTableau(::SDIRK_Midpoint_1_2, ::Type{T}=Float64) where {T}
  SDIRK12(T)
end

############
# 2 stages #
############
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

"""
    DIRK_CrankNicolson_2_2
"""
struct DIRK_CrankNicolson_2_2 <: TableauName end

function ButcherTableau(::DIRK_CrankNicolson_2_2, ::Type{T}=Float64) where {T}
  DIRK22(0, 1, 1 // 2, T)
end

"""
    SDIRK_QinZhang_2_2
"""
struct SDIRK_QinZhang_2_2 <: TableauName end

function ButcherTableau(::SDIRK_QinZhang_2_2, ::Type{T}=Float64) where {T}
  DIRK22(1 // 4, 3 // 4, 1 // 4, T)
end

"""
    DIRK_LobattoIIIA_2_2
"""
struct DIRK_LobattoIIIA_2_2 <: TableauName end

function ButcherTableau(::DIRK_LobattoIIIA_2_2, ::Type{T}=Float64) where {T}
  tableau = DIRK22(0, 1, 1 // 2, T)
  emb_weights = T[1, 0]
  emb_order = 1
  EmbeddedTableau(tableau, emb_weights, emb_order)
end

"""
    DIRK_RadauI_2_3
"""
struct DIRK_RadauI_2_3 <: TableauName end

function ButcherTableau(::DIRK_RadauI_2_3, ::Type{T}=Float64) where {T<:Real}
  a = 1 // 3
  b = 1 // 4
  c = 3 // 4
  matrix = T[
    0 0
    a a
  ]
  weights = T[b, c]
  order = 3
  GenericTableau(matrix, weights, order)
end

"""
    DIRK_RadauII_2_3
"""
struct DIRK_RadauII_2_3 <: TableauName end

function ButcherTableau(::DIRK_RadauII_2_3, ::Type{T}=Float64) where {T<:Real}
  a = 1 // 3
  b = 3 // 4
  c = 1 // 4
  matrix = T[
    a 0
    1 0
  ]
  weights = T[b, c]
  order = 3
  GenericTableau(matrix, weights, order)
end

"""
    SDIRK_LobattoIIIC_2_2
"""
struct SDIRK_LobattoIIIC_2_2 <: TableauName end

function ButcherTableau(::SDIRK_LobattoIIIC_2_2, ::Type{T}=Float64) where {T}
  DIRK22(0, 1, 0, T)
end

"""
    SDIRK_2_2
"""
struct SDIRK_2_2 <: TableauName end

function ButcherTableau(::SDIRK_2_2, ::Type{T}=Float64) where {T}
  DIRK22(1, 0, 1, T)
end

"""
    SDIRK_SSP_2_3
    SDIRK_Crouzeix_2_3
"""
struct SDIRK_SSP_2_3 <: TableauName end

function ButcherTableau(::SDIRK_SSP_2_3, ::Type{T}=Float64) where {T}
  DIRK23(-1, T)
end

struct SDIRK_Crouzeix_2_3 <: TableauName end

function ButcherTableau(::SDIRK_Crouzeix_2_3, ::Type{T}=Float64) where {T}
  ButcherTableau(SDIRK_SSP_2_3(), T)
end

############
# 3 stages #
############
"""
    SDIRK_3_2
"""
struct SDIRK_3_2 <: TableauName end

function ButcherTableau(::SDIRK_3_2, ::Type{T}=Float64) where {T<:Real}
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
    DIRK_TRBDF_3_2
"""
struct DIRK_TRBDF_3_2 <: TableauName end

function ButcherTableau(::DIRK_TRBDF_3_2, ::Type{T}=Float64) where {T<:Real}
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
    DIRK_TRX_3_2
"""
struct DIRK_TRX_3_2 <: TableauName end

function ButcherTableau(::DIRK_TRX_3_2, ::Type{T}=Float64) where {T<:Real}
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

"""
    SDIRK_3_3
"""
struct SDIRK_3_3 <: TableauName end

function ButcherTableau(::SDIRK_3_3, ::Type{T}=Float64) where {T}
  α = 2 * cospi(1 // 18) / sqrt(3)
  a = (1 - α) / 2
  b = -3 * α^2 + 4 * α - 1 // 4
  c = 3 * α^2 - 5 * α + 5 // 4
  matrix = T[
    α 0 0
    a α 0
    b c α
  ]
  weights = T[b, c, α]
  order = 3
  GenericTableau(matrix, weights, order)
end

"""
    SDIRK_Crouzeix_3_4
"""
struct SDIRK_Crouzeix_3_4 <: TableauName end

function ButcherTableau(::SDIRK_Crouzeix_3_4, ::Type{T}=Float64) where {T}
  α = 2 * cospi(1 // 18) / sqrt(3)
  a = (1 + α) / 2
  b = -α / 2
  c = 1 + α
  d = -(1 + 2 * α)
  e = 1 // 6 / α^2
  f = 1 - 1 // 3 / α^2
  matrix = T[
    a 0 0
    b a 0
    c d a
  ]
  weights = T[e, f, e]
  order = 4
  GenericTableau(matrix, weights, order)
end

"""
    SDIRK_Norsett_3_4
"""
struct SDIRK_Norsett_3_4 <: TableauName end

function ButcherTableau(::SDIRK_Norsett_3_4, ::Type{T}=Float64) where {T}
  # One of the three roots of x^3 - 3*x^2 + x/2 - 1/24
  # The largest one brings most stability
  α = 1.0685790213016289
  b = 1 // 2 - α
  c = 2 * α
  d = 1 - 4 * α
  e = 1 // 6 / (1 - 2 * α)^2
  f = 1 - 1 // 3 / (1 - 2 * α)^2
  matrix = T[
    α 0 0
    b α 0
    c d α
  ]
  weights = T[e, f, e]
  order = 4
  GenericTableau(matrix, weights, order)
end

"""
    DIRK_LobattoIIIC_3_4
"""
struct DIRK_LobattoIIIC_3_4 <: TableauName end

function ButcherTableau(::DIRK_LobattoIIIC_3_4, ::Type{T}=Float64) where {T}
  a = 1 // 4
  b = 1 // 6
  c = 2 // 3
  matrix = T[
    0 0 0
    a a 0
    0 1 0
  ]
  weights = T[b, c, b]
  order = 4
  GenericTableau(matrix, weights, order)
end

############
# 4 stages #
############
"""
    SDIRK_4_3
"""
struct SDIRK_4_3 <: TableauName end

function ButcherTableau(::SDIRK_4_3, ::Type{T}=Float64) where {T}
  a = 1 // 2
  b = 1 // 6
  c = -1 // 2
  d = 3 // 2
  e = -3 // 2
  matrix = T[
    a 0 0 0
    b a 0 0
    c a a 0
    d e a a
  ]
  weights = T[d, e, a, a]
  order = 4
  GenericTableau(matrix, weights, order)
end
