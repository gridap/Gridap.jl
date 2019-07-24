module GridGraphs

using Test
using Gridap
using Gridap.Helpers

export GridGraph
export FullGridGraph
export cell_nfaces
export nface_cells
export connections
export test_grid_graph
export test_full_grid_graph
import Base: ndims

# Interfaces

"""
Abstract type that provides aggregation info of cell to nfaces and nface to cells
"""
abstract type GridGraph end

function ndims(::GridGraph)
  @abstractmethod
end

function cell_nfaces(g::GridGraph,n::Integer)::IndexCellVector
  @abstractmethod
end

function nface_cells(g::GridGraph,n::Integer)::IndexCellVector
  @abstractmethod
end

function connections(g::GridGraph,from::Integer,to::Integer)
  nd = ndims(g)
  if from == nd
    return cell_nfaces(g,to)
  elseif to == nd
    return nface_cells(g,from)
  else
    @unreachable
  end
end

"""
Abstract type that provides aggregation info of nfaces to mfaces for all n,m
"""
abstract type FullGridGraph <: GridGraph end

function connections(g::FullGridGraph,from::Integer,to::Integer)::IndexCellVector
  @abstractmethod
end

function cell_nfaces(g::FullGridGraph,n::Integer)
  from = ndims(g)
  connections(g,from,n)
end

function nface_cells(g::FullGridGraph,n::Integer)
  to = ndims(g)
  connections(g,n,to)
end

# Testers

function test_grid_graph(g::GridGraph,nd::Integer)
  _nd = ndims(g)
  @test isa(_nd,Integer)
  @test _nd == nd
  for d in 0:nd
    ca = cell_nfaces(g,d)
    cb = connections(g,nd,d)
    @test isa(ca,IndexCellVector{<:Integer})
    @test isa(cb,IndexCellVector{<:Integer})
    @test ca == cb
    ca = nface_cells(g,d)
    cb = connections(g,d,nd)
    @test isa(ca,IndexCellVector{<:Integer})
    @test isa(cb,IndexCellVector{<:Integer})
    @test ca == cb
  end
end

function test_full_grid_graph(g::FullGridGraph,nd::Integer)
  test_grid_graph(g,nd)
  for i in 0:nd
    for j in 0:nd
      ca = connections(g,i,j)
      @test isa(ca,IndexCellVector{<:Integer})
    end
  end
end

# Concrete implementations

struct GridGraphFromData <: GridGraph
  primal::Vector{IndexCellVector}
  dual::Vector{IndexCellVector}
end

function ndims(g::GridGraphFromData)
  n = length(g.primal)
  m = length(g.dual)
  @assert n == m
  n-1
end

function cell_nfaces(g::GridGraphFromData,n::Integer)
  g.primal[n+1]
end

function nface_cells(g::GridGraphFromData,n::Integer)
  g.dual[n+1]
end

struct FullGridGraphFromData <: FullGridGraph
  data::Array{IndexCellArray,2}
end

function ndims(g::FullGridGraphFromData)
  n,m = size(g.data)
  @assert n == m
  n-1
end

function connections(g::FullGridGraphFromData,from::Integer,to::Integer)
  g.data[from+1,to+1]
end

end # module
