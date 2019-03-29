module Quadratures

using Numa # @santiagobadia : To be eliminated
using Numa.Helpers
using Numa.FieldValues
using QuadGK

export Quadrature, coordinates, weights
export TensorProductQuadrature, TensorProductQuadratureOld

# Abstract types and interfaces

"""
Abstract type representing a quadrature rule on a Polytope in a space of D dimensions
"""
abstract type Quadrature{D} end

coordinates(::Quadrature{D} where D)::Array{Point{D},1} = @abstractmethod

weights(::Quadrature)::Array{Float64,1} = @abstractmethod

#Concrete structs

"""
Tensor product quadrature rule (nodes and weights) on a hyper cube [-1,1]^D
"""
struct TensorProductQuadrature{D} <: Quadrature{D}
  coords::Array{Point{D},1}
  weights::Array{Float64,1}
end

# Methods

include("QuadraturesMethods.jl")

end # module Quadratures
