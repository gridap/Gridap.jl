module Grids

using Test
using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery

export Grid
export points
export cells
export celltypes
export cellorders
export test_grid
export celldim
export pointdim
import Gridap: ncells
export npoints
export CellPolytopes
import Gridap: Triangulation
import Gridap: CellPoints
import Gridap: CellRefFEs
import Gridap: CellBasis

# Interfaces

"""
Abstract type representing a FE mesh a.k.a. grid (i.e., there is a 1 to 1
correspondence between a grid and a conforming Lagrangian FE space with no
Dirichlet boundary conditions).
D is the dimension of the coordinates and Z is the dimension of the cells.
"""
abstract type Grid{D,Z} end

function points(::Grid{D})::IndexCellValue{<:Point{D}} where D
  @abstractmethod
end

cells(::Grid)::IndexCellArray{Int,1} = @abstractmethod

function celltypes(::Grid{D,Z})::CellValue{NTuple{Z,Int}} where {D,Z}
  @abstractmethod
end

cellorders(::Grid)::CellValue{Int} = @abstractmethod

celldim(::Grid{D,Z}) where {D,Z} = Z

pointdim(::Grid{D,Z}) where {D,Z} = D

ncells(g::Grid) = length(celltypes(g))

npoints(g::Grid) = length(points(g))

CellPolytopes(g::Grid) = _cell_polytopes(celltypes(g))

Triangulation(grid::Grid) = TriangulationFromGrid(grid)

# Testers

function test_grid(grid::Grid{D,Z},np,nc) where {D,Z}
  @test isa(points(grid),IndexCellValue{<:Point{D}})
  @test isa(cells(grid),IndexCellVector{<:Integer})
  @test isa(celltypes(grid),CellValue{NTuple{Z,Int}})
  @test isa(cellorders(grid),IndexCellValue{Int})
  @test np == length(points(grid))
  @test nc == length(cells(grid))
  @test nc == length(celltypes(grid))
  @test nc == length(cellorders(grid))
  trian = Triangulation(grid)
  test_triangulation(trian)
end

# Pretty printing

import Base: show

function show(io::IO,self::Grid{D,Z}) where {D,Z}
  print(io,"$(nameof(typeof(self))) object")
end

function show(io::IO,::MIME"text/plain",grid::Grid)
  show(io,grid)
  print(io,":")
  print(io,"\n pointdim: $(pointdim(grid))")
  print(io,"\n celldim: $(celldim(grid))")
  print(io,"\n npoints: $(npoints(grid))")
  print(io,"\n ncells: $(ncells(grid))")
end

# Concrete implementation

struct GridFromData{D,Z} <: Grid{D,Z}
  _points::IndexCellValue{<:Point{D}}
  _cells::IndexCellVector
  _celltypes::IndexCellValue{NTuple{Z,Int}}
  _cellorders::IndexCellValue{Int}
end

points(g::GridFromData) = g._points

cells(g::GridFromData) = g._cells

celltypes(g::GridFromData) = g._celltypes

cellorders(g::GridFromData) = g._cellorders

# Helpers

struct TriangulationFromGrid{D,Z,G<:Grid{D,Z}} <: Triangulation{Z,D}
  grid::G
end

function CellPoints(self::TriangulationFromGrid)
  CellVectorFromLocalToGlobal(cells(self.grid),points(self.grid))
end

function CellRefFEs(trian::TriangulationFromGrid)
  grid = trian.grid
  _build_cell_ref_fes(celltypes(grid),cellorders(grid))
end

function CellBasis(trian::TriangulationFromGrid)
  reffes = CellRefFEs(trian)
  _setup_cell_basis(reffes)
end

function _build_cell_ref_fes(codes,orders)
  @notimplemented
end

function _build_cell_ref_fes(
  codes::ConstantCellValue, orders::ConstantCellValue)
  code = codes.value
  D = length(code)
  order = orders.value
  polytope = Polytope(code)
  _orders = fill(order,D)
  reffe = LagrangianRefFE(Float64,polytope,_orders)
  ConstantCellValue(reffe,length(codes))
end

function _setup_cell_basis(reffes)
  @notimplemented
end

function _setup_cell_basis(reffes::ConstantCellValue)
  reffe = reffes.value
  basis = shfbasis(reffe)
  ConstantCellMap(basis,reffes.length)
end

function _cell_polytopes(ct)
  @notimplemented
end

function _cell_polytopes(ct::ConstantCellValue)
  extr = ct.value
  l = ct.length
  p = Polytope(extr)
  ConstantCellValue(p,l)
end

end # module
