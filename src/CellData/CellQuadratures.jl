"""
    struct CellQuadrature <: GridapType
      # + private fields
    end
"""
struct CellQuadrature <: GridapType
  array
  coords
  weights
end

"""
    CellQuadrature(array::AbstractArray{<:Quadrature})
"""
function CellQuadrature(array::AbstractArray{<:Quadrature})
  coords = _get_coordinates(array)
  weights = _get_weights(array)
  CellQuadrature(array,coords,weights)
end

"""
    CellQuadrature(coords::AbstractArray,weights::AbstractArray)
"""
function CellQuadrature(coords,weights)
  array = apply((a,b)->GenericQuadrature(a,b),coords,weights)
  CellQuadrature(array,coords,weights)
end

"""
    CellQuadrature(degree,polytopes::Vector{<:Polytope}, cell_types::AbstractVector)
"""
function CellQuadrature(degree,polytopes::Vector{<:Polytope}, cell_types::AbstractVector)
  f = (p) -> Quadrature(p,degree)
  quads = map(f,polytopes)
  array = CompressedArray(quads,cell_types)
  CellQuadrature(array)
end

function CellQuadrature(degree,polytopes::Vector{<:Polytope}, cell_types::Fill)
  ctype = cell_types.value
  p = polytopes[ctype]
  quad = Quadrature(p,degree)
  array = Fill(quad,length(cell_types))
  CellQuadrature(array)
end

"""
    get_array(quad::CellQuadrature)
"""
get_array(quad::CellQuadrature) = quad.array

"""
    get_coordinates(q::CellQuadrature)
"""
function get_coordinates(quad::CellQuadrature)
  quad.coords
end

"""
    get_weights(q::CellQuadrature)
"""
function get_weights(quad::CellQuadrature)
  quad.weights
end

function _get_coordinates(q::AbstractArray{<:Quadrature})
  @notimplemented "Not implemented, since we dont need it"
end

function _get_coordinates(q::CompressedArray{<:Quadrature})
  coords = map(get_coordinates,q.values)
  CompressedArray(coords,q.ptrs)
end

function _get_coordinates(q::Fill{<:Quadrature})
  coords = get_coordinates(q.value)
  Fill(coords,length(q))
end

function _get_coordinates(q::AppendedArray)
  a = _get_coordinates(q.a)
  b = _get_coordinates(q.b)
  lazy_append(a,b)
end

function _get_weights(q::AbstractArray{<:Quadrature})
  @notimplemented "Not implemented, since we dont need it"
end

function _get_weights(q::CompressedArray{<:Quadrature})
  w = map(get_weights,q.values)
  CompressedArray(w,q.ptrs)
end

function _get_weights(q::Fill{<:Quadrature})
  w = get_weights(q.value)
  Fill(w,length(q))
end

function _get_weights(q::AppendedArray)
  a = _get_weights(q.a)
  b = _get_weights(q.b)
  lazy_append(a,b)
end

"""
    integrate(cell_field,ϕ::CellField,quad::CellQuadrature)
"""
function integrate(f,ϕ::CellField,quad::CellQuadrature)
  integrate(f,ϕ,get_coordinates(quad),get_weights(quad))
end

function integrate(cell_field,ϕ::CellField,q,w)
  f = convert_to_cell_field(cell_field,length(ϕ))
  j = get_array(∇(ϕ))
  @assert length(f) == length(ϕ) "Are you using the right triangulation to integrate?"
  @assert length(f) == length(w) "Are you using the right quadrature to integrate?"
  integrate(get_array(f∘ϕ),q,w,j)
end

function lazy_append(quad1::CellQuadrature,quad2::CellQuadrature)
  array = lazy_append(quad1.array,quad2.array)
  coords = lazy_append(quad1.coords,quad2.coords)
  weights = lazy_append(quad1.weights,quad2.weights)
  CellQuadrature(array,coords,weights)
end

# Some syntactic sugar

struct MappedCellValues{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
  physvals::A
  refvals::AbstractArray
  ϕ::CellField
end

Base.size(a::MappedCellValues) = size(a.physvals)
Base.getindex(a::MappedCellValues,i...) = a.physvals[i...]
Base.IndexStyle(::Type{MappedCellValues{T,N,A}}) where {T,N,A} = IndexStyle(A)
Arrays.get_array(a::MappedCellValues) = a.physvals
Arrays.array_cache(a::MappedCellValues) = array_cache(a.physvals)
Arrays.getindex!(cache,a,i...) = getindex!(cache,a.physvals,i...)

function evaluate(f::CellField,x::MappedCellValues)
  q = x.refvals
  ϕ = x.ϕ
  evaluate(f∘ϕ,q)
end

function evaluate(ϕ::CellField,quad::CellQuadrature)
  q = get_coordinates(quad)
  wq = get_weights(quad)
  x = MappedCellValues(ϕ(q),q,ϕ)
  wx = MappedCellValues( apply(bcast(*),wq, det(∇(ϕ))(q)), wq, ϕ)
  CellQuadrature(x,wx)
end

function integrate(f,quad::CellQuadrature)
  x = get_coordinates(quad)
  wx = get_weights(quad)
  _integrate(f,x,wx)
end

function  _integrate(f,x,wx)
  @notimplemented
end

function  _integrate(f,x::MappedCellValues,wx::MappedCellValues)
  @assert x.ϕ == wx.ϕ
  ϕ = x.ϕ
  q = x.refvals
  wq = wx.refvals
  integrate(f,ϕ,q,wq)
end

struct Integral
  integrand
end

const ∫ = Integral

function Base.:*(a::Integral,b::CellQuadrature)
  integrate(a.integrand,b)
end

