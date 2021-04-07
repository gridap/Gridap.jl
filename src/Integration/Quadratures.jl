
"""
    abstract type Quadrature{D,T} <: GridapType end

-[`get_coordinates(q::Quadrature)`](@ref)
-[`get_weights(q::Quadrature)`](@ref)
-[`test_quadrature`](@ref)
"""
abstract type Quadrature{D,T} <: GridapType end

# Abstract API

"""
    get_coordinates(q::Quadrature)
"""
function get_coordinates(q::Quadrature)
  @abstractmethod
end

"""
    get_weights(q::Quadrature)
"""
function get_weights(q::Quadrature)
  @abstractmethod
end

# Default API

"""
    num_points(q::Quadrature)
"""
num_points(q::Quadrature) = length(get_weights(q))

"""
    num_point_dims(::Quadrature{D}) where D
    num_point_dims(::Type{<:Quadrature{D}}) where D
"""
num_point_dims(::Quadrature{D}) where D = D

num_point_dims(::Type{<:Quadrature{D}}) where D = D

"""
    num_dims(::Quadrature{D}) where D where D
    num_dims(::Type{<:Quadrature{D}}) where D
"""
num_dims(::Quadrature{D}) where D = D

num_dims(::Type{<:Quadrature{D}}) where D = D

# Tester

"""
    test_quadrature(q::Quadrature{D,T}) where {D,T}
"""
function test_quadrature(q::Quadrature{D,T}) where {D,T}
  x = get_coordinates(q)
  w = get_weights(q)
  @test isa(x,Vector{Point{D,T}})
  @test isa(w,Vector{T})
  @test length(x) == num_points(q)
  @test length(w) == num_points(q)
  @test D == num_dims(q)
  @test D == num_point_dims(q)
  @test D == num_dims(typeof(q))
  @test D == num_point_dims(typeof(q))
end

# Generic concrete implementation

"""
    struct GenericQuadrature{D,T} <: Quadrature{D,T}
      coordinates::Vector{Point{D,T}}
      weights::Vector{T}
    end
"""
struct GenericQuadrature{D,T} <: Quadrature{D,T}
  coordinates::Vector{Point{D,T}}
  weights::Vector{T}
end

get_coordinates(q::GenericQuadrature) = q.coordinates

get_weights(q::GenericQuadrature) = q.weights

function GenericQuadrature(a::Quadrature)
  GenericQuadrature(get_coordinates(q),get_weights(a))
end

function GenericQuadrature(a::GenericQuadrature)
  a
end
