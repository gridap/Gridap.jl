
"""
    Quadrature(polytope::Polytope{D},degree) where D
"""
function Quadrature(p::Polytope{D},degree) where D
  if is_n_cube(p)
    q = TensorProductQuadrature{D}(degree)
    quad = GenericQuadrature(q.coordinates,q.weights)
  elseif is_simplex(p)
    q = DuffyQuadrature{D}(degree)
    quad = GenericQuadrature(q.coordinates,q.weights)
  else
    @notimplemented "Quadratures only implemented for n-cubes and simplices"
  end
  quad
end

"""
    CellQuadrature(polytopes::Vector{<:Polytope}, cell_types::AbstractVector)
"""
function CellQuadrature(degree,polytopes::Vector{<:Polytope}, cell_types::AbstractVector)
  f = (p) -> Quadrature(p,degree)
  quads = map(f,polytopes)
  CompressedArray(quads,cell_types)
end

function CellQuadrature(degree,polytopes::Vector{<:Polytope}, cell_types::Fill)
  ctype = cell_types.value
  p = polytopes[ctype]
  quad = Quadrature(p,degree)
  Fill(quad,length(cell_types))
end

"""
    CellQuadrature(trian::Triangulation, degree)
"""
function CellQuadrature(trian::Triangulation, degree)
  polytopes = map(get_polytope,get_reffes(trian))
  cell_type = get_cell_type(trian)
  CellQuadrature(degree,polytopes,cell_type)
end

"""
    get_coordinates(q::AbstractArray{<:Quadrature})
"""
function get_coordinates(q::AbstractArray{<:Quadrature})
  @notimplemented "Not implemented, since we dont need it"
end

function get_coordinates(q::CompressedArray{<:Quadrature})
  coords = map(get_coordinates,q.values)
  CompressedArray(coords,q.ptrs)
end

function get_coordinates(q::Fill{<:Quadrature})
  coords = get_coordinates(q.value)
  Fill(coords,length(q))
end

"""
    get_weights(q::AbstractArray{<:Quadrature})
"""
function get_weights(q::AbstractArray{<:Quadrature})
  @notimplemented "Not implemented, since we dont need it"
end

function get_weights(q::CompressedArray{<:Quadrature})
  w = map(get_weights,q.values)
  CompressedArray(w,q.ptrs)
end

function get_weights(q::Fill{<:Quadrature})
  w = get_weights(q.value)
  Fill(w,length(q))
end

"""
the `cell_field` is aligned with the cells in `trian`
"""
function integrate(cell_field,trian::Triangulation,quad::AbstractArray{<:Quadrature})
  cell_map = get_cell_map(trian)
  q = get_coordinates(quad)
  w = get_weights(quad)
  j = gradient(cell_map)
  f = convert_to_cell_field(cell_field,cell_map)
  integrate(get_array(f),q,w,j)
end

