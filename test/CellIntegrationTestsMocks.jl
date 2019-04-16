
# TODO: the following is only a dummy implementation in order to
# test the interfaces at this early stage.
# It will be removed in the future and replaced by other more general
# and elegant implementations

# include("PolynomialsTestsMocks.jl")

using Numa.CellValues: IndexCellArray
import Numa.CellValues: cellsize
using Numa.CellIntegration
import Numa.CellIntegration: cellcoordinates, cellbasis
using Numa.Polytopes
using Numa.RefFEs
using Numa.FieldValues
import Numa.Geometry: celltypes

struct DummyCellCoordinates2D <: IndexCellArray{Point{2},1,Array{Point{2},1},1}
  x::Array{Point{2},2}
  c::Array{Point{2},1}
end

function DummyCellCoordinates2D(;partition::Tuple{Int,Int})
  ncells = prod(partition)
  x = Array{Point{2},2}(undef,(4,ncells))
  pc = MPoint{2}(zeros(2))
  pn = MPoint{2}(zeros(2))
  hx = 1.0/partition[1]
  hy = 1.0/partition[2]
  cell = 1
  for cy in 1:partition[2]
    pc[2] = (cy-1)*hy
    for cx in 1:partition[1]
      pc[1] = (cx-1)*hx
      node = 1
      for ny in 1:2
        pn[2] = pc[2] + (ny-1)*hy
        for nx in 1:2
          pn[1] = pc[1] + (nx-1)*hx
          x[node,cell] = pn
          node += 1
        end
      end
      cell += 1
    end
  end
  c = Array{Point{2},1}(undef,(4,))
  DummyCellCoordinates2D(x,c)
end

function Base.getindex(self::DummyCellCoordinates2D,cell::Int)
  @inbounds for i in 1:4
    self.c[i] = self.x[i,cell]
  end
  self.c
end

Base.size(self::DummyCellCoordinates2D) = (size(self.x,2),)

IndexStyle(::Type{DummyCellCoordinates2D}) = IndexLinear()

cellsize(self::DummyCellCoordinates2D) = (4,)

struct DummyIntegrationMesh2D <: Triangulation{2,2}
  cellcoords::DummyCellCoordinates2D
  cellbasis::CellBasisFromSingleInterpolation{2,Float64}
  celltypes::ConstantCellValue{NTuple{2,Int}}
end

function DummyIntegrationMesh2D(;partition::Tuple{Int,Int})
  # basis = ShapeFunctionsScalarQua4()
  polytope = Polytope(Polytopes.PointInt{2}(1,1))
  reffe = LagrangianRefFE{2,ScalarValue}(polytope,[1,1])
  basis = reffe.shfbasis
  cellbasis = CellBasisFromSingleInterpolation(basis)
  cellcoords = DummyCellCoordinates2D(partition=partition)
  celltypes = ConstantCellValue((HEX_AXIS,HEX_AXIS),length(cellcoords))
  DummyIntegrationMesh2D(cellcoords,cellbasis,celltypes)
end

cellcoordinates(self::DummyIntegrationMesh2D) = self.cellcoords

cellbasis(self::DummyIntegrationMesh2D) = self.cellbasis

celltypes(self::DummyIntegrationMesh2D) = self.celltypes
