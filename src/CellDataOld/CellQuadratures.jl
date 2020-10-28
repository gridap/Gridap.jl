"""
    struct CellQuadrature <: GridapType
      array
    end
"""
struct CellQuadrature <: GridapType
  array
  coords
  weights
  @doc """
      CellQuadrature(array::AbstractArray{<:Quadrature})
  """
  function CellQuadrature(array::AbstractArray{<:Quadrature})
    coords = _get_coordinates(array)
    weights = _get_weights(array)
    new(array,coords,weights)
  end
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
    integrate(cell_field,cell_map::AbstractArray{<:Field},quad::CellQuadrature)
"""
function integrate(cell_field,cell_map::AbstractArray{<:Field},quad::CellQuadrature)
  q = get_coordinates(quad)
  w = get_weights(quad)
  j = gradient(cell_map)
  _f = convert_to_cell_field(cell_field,cell_map)
  f = to_ref_space(_f)
  @assert length(f) == length(cell_map) "Are you using the right triangulation to integrate?"
  @assert length(f) == length(w) "Are you using the right quadrature to integrate?"
  integrate(get_array(f),q,w,j)
end

function lazy_append(quad1::CellQuadrature,quad2::CellQuadrature)
  array = lazy_append(quad1.array,quad2.array)
  CellQuadrature(array)
end

