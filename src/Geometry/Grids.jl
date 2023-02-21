
"""
    abstract type Grid{Dc,Dp}

Abstract type that represents mesh of a domain of parametric dimension `Dc` and
physical dimension `Dp`.

The interface of `Grid` is defined by overloading:

- [`get_node_coordinates(trian::Grid)`](@ref)
- [`get_cell_node_ids(trian::Grid)`](@ref)
- [`get_reffes(trian::Grid)`](@ref)
- [`get_cell_type(trian::Grid)`](@ref)

The `Grid`  interface has the following traits

- [`OrientationStyle(::Type{<:Grid})`](@ref)
- [`RegularityStyle(::Type{<:Grid})`](@ref)

The interface of `Grid` is tested with

- [`test_grid`](@ref)

"""
abstract type Grid{Dc,Dp} <: GridapType end

# Traits

abstract type OrientationStyle end
struct Oriented <: OrientationStyle end
struct NonOriented <: OrientationStyle end

"""
    OrientationStyle(::Type{<:Grid})
    OrientationStyle(::Grid)

`Oriented()` if has oriented faces, `NonOriented()` otherwise (default).
"""
OrientationStyle(a::Grid) = OrientationStyle(typeof(a))
OrientationStyle(::Type{<:Grid}) = NonOriented()

abstract type RegularityStyle end
struct Regular <: RegularityStyle end
struct Irregular <: RegularityStyle end

"""
    RegularityStyle(::Type{<:Grid})
    RegularityStyle(::Grid)

`Regular()` if no hanging-nodes default), `Irregular()` otherwise.
"""
RegularityStyle(a::Grid) = RegularityStyle(typeof(a))
RegularityStyle(::Type{<:Grid}) = Regular()

# Interface

"""
    get_node_coordinates(trian::Grid) -> AbstractArray{<:Point{Dp}}
"""
function get_node_coordinates(trian::Grid)
  @abstractmethod
end

"""
    get_cell_node_ids(trian::Grid)
"""
function get_cell_node_ids(trian::Grid)
  @abstractmethod
end

"""
    get_reffes(trian::Grid) -> Vector{LagrangianRefFE}
"""
function get_reffes(trian::Grid)
  @abstractmethod
end

"""
    get_cell_type(trian::Grid) -> AbstractVector{<:Integer}
"""
function get_cell_type(trian::Grid)
  @abstractmethod
end

"""
    get_facet_normal(trian::Grid)
"""
function get_facet_normal(trian::Grid)
  Dp = num_point_dims(trian)
  Dc = num_cell_dims(trian)
  if Dp == Dc + 1
    @abstractmethod
  else
    @unreachable "get_facet_normal does not make sense for this grid"
  end
end

"""
    test_grid(trian::Grid)
"""
function test_grid(trian::Grid{Dc,Dp}) where {Dc,Dp}
  @test num_cell_dims(trian) == Dc
  @test num_point_dims(trian) == Dp
  @test num_cell_dims(typeof(trian)) == Dc
  @test num_point_dims(typeof(trian)) == Dp
  cell_coords = get_cell_coordinates(trian)
  @test isa(cell_coords,AbstractArray{<:AbstractVector{<:Point}})
  reffes = get_reffes(trian)
  @test isa(reffes,AbstractVector{<:LagrangianRefFE{Dc}})
  cell_types = get_cell_type(trian)
  @test isa(cell_types,AbstractArray{<:Integer})
  ncells = num_cells(trian)
  @test ncells == length(cell_coords)
  @test ncells == length(cell_types)
  nodes_coords = get_node_coordinates(trian)
  @test isa(nodes_coords,AbstractArray{<:Point})
  cell_node_ids = get_cell_node_ids(trian)
  @test isa(cell_node_ids,AbstractArray{<:AbstractArray{<:Integer}})
  @test num_nodes(trian) == length(nodes_coords)
  @test isa(is_oriented(trian),Bool)
  @test isa(is_regular(trian),Bool)
  @test OrientationStyle(trian) in (Oriented(), NonOriented())
  @test RegularityStyle(trian) in (Regular(), Irregular())
  @test is_oriented(trian) == (OrientationStyle(trian) == Oriented())
  @test is_regular(trian) == (RegularityStyle(trian) == Regular())
end

# Some API

"""
    num_cells(trian::Grid) -> Int
"""
num_cells(trian::Grid) = length(get_cell_type(trian))

"""
    num_nodes(trian::Grid) -> Int
"""
num_nodes(trian::Grid) = length(get_node_coordinates(trian))

"""
    num_cell_dims(::Grid) -> Int
    num_cell_dims(::Type{<:Grid}) -> Int
"""
num_cell_dims(::Grid{Dc,Dp}) where {Dc,Dp} = Dc
num_cell_dims(::Type{<:Grid{Dc,Dp}}) where {Dc,Dp} = Dc

"""
    num_point_dims(::Grid) -> Int
    num_point_dims(::Type{<:Grid}) -> Int
"""
num_point_dims(::Grid{Dc,Dp}) where {Dc,Dp} = Dp
num_point_dims(::Type{<:Grid{Dc,Dp}}) where {Dc,Dp} = Dp

"""
    num_dims(::Grid) -> Int
    num_dims(::Type{<:Grid}) -> Int

Equivalent to `num_cell_dims`.
"""
num_dims(g::Grid{Dc}) where Dc = Dc
num_dims(::Type{<:Grid{Dc}}) where Dc = Dc

"""
    is_first_order(trian::Grid) -> Bool
"""
function is_first_order(trian::Grid)
  reffes = get_reffes(trian)
  all(map(is_first_order,reffes))
end

"""
    get_cell_reffe(trian::Grid) -> Vector{<:LagrangianRefFE}

It is not desirable to iterate over the resulting array
for large number of cells if the underlying reference FEs
are of different Julia type.
"""
function get_cell_reffe(trian::Grid)
  type_to_reffe = get_reffes(trian)
  cell_to_type = get_cell_type(trian)
  expand_cell_data(type_to_reffe,cell_to_type)
end

"""
"""
function get_cell_ref_coordinates(trian::Grid)
  type_to_reffe = get_reffes(trian)
  type_to_coords = map(get_node_coordinates,type_to_reffe)
  cell_to_type = get_cell_type(trian)
  expand_cell_data(type_to_coords,cell_to_type)
end

"""
    get_cell_shapefuns(trian::Grid) -> Vector{<:Field}
"""
function get_cell_shapefuns(trian::Grid)
  type_to_reffes = get_reffes(trian)
  cell_to_type = get_cell_type(trian)
  type_to_shapefuns = map(get_shapefuns, type_to_reffes)
  expand_cell_data(type_to_shapefuns,cell_to_type)
end

"""
    get_cell_map(trian::Grid) -> Vector{<:Field}
"""
function get_cell_map(trian::Grid)
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_shapefuns = get_cell_shapefuns(trian)
  lazy_map(linear_combination,cell_to_coords,cell_to_shapefuns)
end

function get_cell_coordinates(trian::Grid)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_node_ids(trian)
  lazy_map(Broadcasting(Reindex(node_to_coords)),cell_to_nodes)
end

function Quadrature(trian::Grid,args...;kwargs...)
  cell_ctype = get_cell_type(trian)
  ctype_polytope = map(get_polytope,get_reffes(trian))
  ctype_quad = map(p->Quadrature(p,args...;kwargs...),ctype_polytope)
  cell_quad = expand_cell_data(ctype_quad,cell_ctype)
end

"""
    is_oriented(::Type{<:Grid}) -> Bool
    is_oriented(a::Grid) -> Bool
"""
is_oriented(a::Grid) = is_oriented(typeof(a))
is_oriented(a::Type{T}) where T<:Grid = OrientationStyle(T) == Oriented()

"""
    is_regular(::Type{<:Grid}) -> Bool
    is_regular(a::Grid) -> Bool
"""
is_regular(a::Grid) = is_regular(typeof(a))
is_regular(a::Type{T}) where T<:Grid = RegularityStyle(T) == Regular()

"""
    Grid(reffe::LagrangianRefFE)
"""
function Grid(reffe::LagrangianRefFE)
  UnstructuredGrid(reffe)
end

"""
    compute_linear_grid(reffe::LagrangianRefFE)
"""
function compute_linear_grid(reffe::LagrangianRefFE)
  p = get_polytope(reffe)
  if get_order(reffe) == 0
    D = num_cell_dims(reffe)
    partition = tfill(1,Val{D}())
  else
    partition = get_orders(reffe)
  end
  compute_reference_grid(p,partition)
end

"""
    compute_reference_grid(p::LagrangianRefFE, nelems::Integer)
"""
function compute_reference_grid(reffe::LagrangianRefFE, nelems::Integer)
  p = get_polytope(reffe)
  compute_reference_grid(p,nelems)
end

"""
    compute_reference_grid(p::LagrangianRefFE, partition::NTuple{D,Integer})
"""
function compute_reference_grid(reffe::LagrangianRefFE{D}, partition::NTuple{D,Integer}) where D
  p = get_polytope(reffe)
  compute_reference_grid(p,partition)
end

"""
    compute_reference_grid(p::Polytope,nelems)
"""
function compute_reference_grid(p::Polytope,nelems)
  if is_n_cube(p)
    _compute_linear_grid_from_n_cube(p,nelems)
  else
    _compute_linear_grid_from_simplex(p,nelems)
  end
end

function _compute_linear_grid_from_n_cube(
  p::Polytope{D},
  partition::NTuple{D,Integer}) where D

  @assert is_n_cube(p)
  pmin, pmax = get_bounding_box(p)
  desc = CartesianDescriptor(pmin,pmax,partition)
  CartesianGrid(desc)
end

function _compute_linear_grid_from_n_cube(
  p::Polytope{D},
  nelems::Integer) where D

  partition = tfill(nelems,Val{D}())
  _compute_linear_grid_from_n_cube(p,partition)
end

function _compute_linear_grid_from_simplex(
  p::Polytope{D},
  partition::NTuple{D,Integer}) where D

  @assert minimum(partition) == maximum(partition) >= 1
  nelems = partition[1]
  _compute_linear_grid_from_simplex(p,nelems)
end

function _compute_linear_grid_from_simplex(p::Polytope{1},nelems::Integer)
  @assert is_simplex(p)
  _compute_linear_grid_coords_from_n_cube(p,nelems)
end

function _compute_linear_grid_from_simplex(p::Polytope,nelems::Integer)
  @assert is_simplex(p)
  X,T = _compute_linear_grid_coords_from_simplex(p,nelems)
  reffes = [ LagrangianRefFE(Float64,p,1) ]
  cell_types = ones(Int8,length(T))
  T = Table(T)
  UnstructuredGrid(X,T,reffes,cell_types)
end

function _compute_linear_grid_coords_from_simplex(p::Polytope{2},n::Integer)
  tri_num(n) = n*(n+1)÷2
  v(n,i,j) = tri_num(n) - tri_num(n-i+1) + j
  @assert is_simplex(p)
  D = 2
  quad_to_tris = ((1,2,3),(2,4,3))
  quad = CartesianIndices( (0:1,0:1) )
  Tp = eltype(get_vertex_coordinates(p))
  n_verts = tri_num(n+1)
  n_cells = tri_num(n)+tri_num(n-1)
  n_verts_x_cell = num_vertices(p)
  X = zeros(Tp,n_verts)
  T = [ zeros(Int,n_verts_x_cell) for i in 1:n_cells ]
  for i in 1:n+1
    for j in 1:n+1-i+1
      vert = v(n+1,i,j)
      X[vert] = Point((i-1)/n,(j-1)/n)
    end
  end
  for i in 1:n
    for j in 1:n-(i-1)
      verts = ntuple( lv-> v(n+1, (i,j).+quad[lv].I ...), Val{2^D}() )
      cell = v(n,i,j)
      T[cell] .= map(Reindex(verts),quad_to_tris[1])
      if (i-1)+(j-1) < n-1
        cell = tri_num(n) + v(n-1,i,j)
        T[cell] .= map(Reindex(verts),quad_to_tris[2])
      end
    end
  end
  X,T
end

function _compute_linear_grid_coords_from_simplex(p::Polytope{3},n::Integer)
  tri_num(n) = n*(n+1)÷2
  tet_num(n) = n*(n+1)*(n+2)÷6
  v(n,i,j) = tri_num(n) - tri_num(n-i+1) + j
  v(n,i,j,k) = tet_num(n) - tet_num(n-i+1) + v(n-i+1,j,k)
  @assert is_simplex(p)
  D = 3
  cube_to_tets = ((1,2,3,5),(2,4,3,6),(3,5,7,6),(2,3,5,6),(3,4,7,6),(4,6,7,8))
  cube = CartesianIndices( (0:1,0:1,0:1) )
  n_core_tets = length(cube_to_tets)-2
  Tp = eltype(get_vertex_coordinates(p))
  n_verts = tet_num(n+1)
  n_cells = tet_num(n)+n_core_tets*tet_num(n-1)+tet_num(n-2)
  n_verts_x_cell = num_vertices(p)
  X = zeros(Tp,n_verts)
  T = [ zeros(Int,n_verts_x_cell) for i in 1:n_cells ]
  for i in 1:n+1
    for j in 1:n+1-(i-1)
      for k in 1:n+1-(i-1)-(j-1)
        vert = v(n+1,i,j,k)
        X[vert] = Point((i-1)/n,(j-1)/n,(k-1)/n)
      end
    end
  end
  for i in 1:n
    for j in 1:n-(i-1)
      for k in 1:n-(i-1)-(j-1)
        verts = ntuple( lv-> v(n+1, (i,j,k).+cube[lv].I ...), Val{2^D}() )
        cell = v(n,i,j,k)
        T[cell] .= map(Reindex(verts),cube_to_tets[1])
        if (i-1)+(j-1)+(k-1) < n-1
          cell = tet_num(n) + (v(n-1,i,j,k)-1)*n_core_tets
          for t in 1:n_core_tets
            T[cell+t] .= map(Reindex(verts),cube_to_tets[t+1])
          end
        end
        if (i-1)+(j-1)+(k-1) < n-2
          cell = tet_num(n) + n_core_tets*tet_num(n-1) + v(n-2,i,j,k)
          T[cell] .= map(Reindex(verts),cube_to_tets[end])
        end
      end
    end
  end
  X,T
end


"""
    Grid(::Type{ReferenceFE{d}},p::Polytope) where d
"""
function Grid(::Type{ReferenceFE{d}},p::Polytope) where d
  UnstructuredGrid(ReferenceFE{d},p)
end

"""
    simplexify(grid::Grid;kwargs...)
"""
function simplexify(grid::Grid;kwargs...)
  simplexify(UnstructuredGrid(grid);kwargs...)
end
