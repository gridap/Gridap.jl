"""
    BernsteinPType{K} <: PolynomialType{K}

Type representing Bernstein polynomials of maximum order `K`
"""
struct BernsteinPType{K} <: PolynomialType{K} end
isHierarchical(::BernsteinPType) = false

struct Bernstein <: PolynomialField end
