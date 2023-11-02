abstract type ButcherTableauType end

struct BE_1_0_1 <: ButcherTableauType end
struct CN_2_0_2 <: ButcherTableauType end
struct SDIRK_2_0_2 <: ButcherTableauType end
struct SDIRK_2_0_3 <: ButcherTableauType end
struct ESDIRK_3_1_2 <: ButcherTableauType end
struct TRBDF2_3_2_3 <: ButcherTableauType end

"""
Butcher tableau
"""
struct ButcherTableau{T <: ButcherTableauType}
  s::Int # stages
  p::Int # embedded order
  q::Int # order
  a::Matrix # A_ij
  b::Vector # b_j
  c::Vector # c_i
  d::Vector # d_j (embedded)
end

# Butcher Tableaus constructors
"""
Backward-Euler

number of stages: 1
embedded method: no
order: 1
"""
function ButcherTableau(::BE_1_0_1)
  s = 1
  p = 0
  q = 1
  a = reshape([1.0],1,1)
  b = [1.0]
  c = [1.0]
  d = [0.0]
  ButcherTableau{BE_1_0_1}(s,p,q,a,b,c,d)
end

"""
Crank-Nicolson (equivalent to trapezoidal rule)

number of stages: 2
embedded method: no
order: 2
"""
function ButcherTableau(type::CN_2_0_2)
  s = 2
  p = 0
  q = 2
  a = [0.0 0.0; 0.5 0.5]
  b = [0.5, 0.5]
  c = [0.0, 1.0]
  d = [0.0, 0.0]
  ButcherTableau{CN_2_0_2}(s,p,q,a,b,c,d)
end

"""
Qin and Zhang's SDIRK

number of stages: 2
embedded method: no
order: 2
"""
function ButcherTableau(type::SDIRK_2_0_2)
  s = 2
  p = 0
  q = 2
  a = [0.25 0.0; 0.5 0.25]
  b = [0.5, 0.5]
  c = [0.25, 0.75]
  d = [0.0, 0.0]
  ButcherTableau{SDIRK_2_0_2}(s,p,q,a,b,c,d)
end

"""
3rd order SDIRK

number of stages: 2
embedded method: no
order: 3
"""
function ButcherTableau(type::SDIRK_2_0_3)
  s = 2
  p = 0
  q = 3
  γ = (3-√(3))/6
  a = [γ 0.0; 1-2γ γ]
  b = [0.5, 0.5]
  c = [γ, 1-γ]
  d = [0.0, 0.0]
  ButcherTableau{SDIRK_2_0_3}(s,p,q,a,b,c,d)
end

function ButcherTableau(type::ESDIRK_3_1_2)
s = 3
p = 1
q = 2
γ = (2-√(2))/2
b₂ = (1 − 2γ)/(4γ)
b̂₂ = γ*(−2 + 7γ − 5(γ^2) + 4(γ^3)) / (2(2γ − 1))
b̂₃ = −2*(γ^2)*(1 − γ + γ^2) / (2γ − 1)
a = [0.0 0.0 0.0; γ γ 0.0; (1 − b₂ − γ) b₂ γ]
b = [(1 − b₂ − γ), b₂, γ]
c = [0.0, 2γ, 1.0]
d = [(1 − b̂₂ − b̂₃), b̂₂, b̂₃]
ButcherTableau{ESDIRK_3_1_2}(s,p,q,a,b,c,d)
end

function ButcherTableau(type::TRBDF2_3_2_3)
  s = 3
  p = 2
  q = 3
  aux = 2.0-√2.0
  a = [0.0 0.0 0.0; aux/2 aux/2 0.0; √2/4 √2/4 aux/2]
  b = [√2/4, √2/4, aux/2]
  c = [0.0, aux, 1.0]
  d = [(1.0-(√2/4))/3, ((3*√2)/4+1.0)/3, aux/6]
  ButcherTableau{TRBDF2_3_2_3}(s,p,q,a,b,c,d)
end

function ButcherTableau(type::Symbol)
  eval(:(ButcherTableau($type())))
end

abstract type IMEXButcherTableauType end

struct IMEX_FE_BE_2_0_1 <: IMEXButcherTableauType end
struct IMEX_Midpoint_2_0_2 <: IMEXButcherTableauType end

"""
Implicit-Explicit Butcher tableaus
"""
struct IMEXButcherTableau{T <: IMEXButcherTableauType}
  s::Int # stages
  p::Int # embedded order
  q::Int # order
  aᵢ::Matrix # A_ij implicit
  aₑ::Matrix # A_ij explicit
  bᵢ::Vector # b_j implicit
  bₑ::Vector # b_j explicit
  c::Vector # c_i
  d::Vector # d_j (embedded)
end

# IMEX Butcher Tableaus constructors
"""
IMEX Forward-Backward-Euler

number of stages: 2
embedded method: no
order: 1
"""
function IMEXButcherTableau(::IMEX_FE_BE_2_0_1)
  s = 2
  p = 0
  q = 1
  aᵢ = [0.0 0.0; 0.0 1.0]
  aₑ = [0.0 0.0; 1.0 0.0]
  bᵢ = [0.0, 1.0]
  bₑ = [0.0, 1.0]
  c = [0.0, 1.0]
  d = [0.0, 0.0]
  IMEXButcherTableau{IMEX_FE_BE_2_0_1}(s,p,q,aᵢ,aₑ,bᵢ,bₑ,c,d)
end

"""
IMEX Midpoint

number of stages: 2
embedded method: no
order: 2
"""
function IMEXButcherTableau(::IMEX_Midpoint_2_0_2)
  s = 2
  p = 0
  q = 2
  aᵢ = [0.0 0.0; 0.0 0.5]
  aₑ = [0.0 0.0; 0.5 0.0]
  bᵢ = [0.0, 1.0]
  bₑ = [0.0, 1.0]
  c = [0.0, 0.5]
  d = [0.0, 0.0]
  IMEXButcherTableau{IMEX_Midpoint_2_0_2}(s,p,q,aᵢ,aₑ,bᵢ,bₑ,c,d)
end


function IMEXButcherTableau(type::Symbol)
  eval(:(IMEXButcherTableau($type())))
end
