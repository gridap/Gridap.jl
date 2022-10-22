
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

"""
    get_name(q::Quadrature)
"""
function get_name(q::Quadrature)
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
  @test isa(x,AbstractVector{Point{D,T}})
  @test isa(w,AbstractVector{T})
  @test length(x) == num_points(q)
  @test length(w) == num_points(q)
  @test D == num_dims(q)
  @test D == num_point_dims(q)
  @test D == num_dims(typeof(q))
  @test D == num_point_dims(typeof(q))
end

# Generic concrete implementation

"""
    struct GenericQuadrature{D,T,C <: AbstractVector{Point{D,T}},W <: AbstractVector{T}} <: Quadrature{D,T}
      coordinates::Vector{Point{D,T}}
      weights::Vector{T}
      name::String
    end
"""
struct GenericQuadrature{D,T,C <: AbstractVector{Point{D,T}},W <: AbstractVector{T}} <: Quadrature{D,T}
  coordinates::C
  weights::W
  name::String
end

function GenericQuadrature(
  coordinates::AbstractVector{Point{D,T}}, weights::AbstractVector{T}) where {D,T}
  name = "Unknown"
  GenericQuadrature(coordinates,weights,name)
end

get_coordinates(q::GenericQuadrature) = q.coordinates

get_weights(q::GenericQuadrature) = q.weights

get_name(q::GenericQuadrature) = q.name

function GenericQuadrature(a::Quadrature)
  GenericQuadrature(get_coordinates(a),get_weights(a),get_name(a))
end

function GenericQuadrature(a::GenericQuadrature)
  a
end

# Quadrature factory

abstract type QuadratureName end

@noinline function Quadrature(p::Polytope,name::QuadratureName,args...;kwargs...)
  @unreachable """\n
  Undefined factory function Quadrature for name $name and the given arguments.
  """
end

Quadrature(name::QuadratureName,args...;kwargs...) = (name, args, kwargs)

"""
    Quadrature(polytope::Polytope{D},degree) where D
"""
function Quadrature(p::Polytope,degree;T::Type{<:AbstractFloat}=Float64)
  if is_n_cube(p)
    quad = Quadrature(p,tensor_product,degree;T=T)
  elseif is_simplex(p)
    D = num_dims(p)
    #if (D==2 && degree in keys(_strang_tri_k2n)) ||
    #   (D==3 && degree in keys(_strang_tet_k2n))
    # For the moment strang only in 3d since
    # there are some 2d strang quadratures that are not accurate
    # as implemented here (to investigate why)
    if (D==3 && degree in keys(_strang_tet_k2n))
      quad = Quadrature(p,strang,degree;T=T)
    else
      quad = Quadrature(p,duffy,degree,T=T)
    end
  else
    @notimplemented "Quadratures only implemented for n-cubes and simplices"
  end
  quad
end
