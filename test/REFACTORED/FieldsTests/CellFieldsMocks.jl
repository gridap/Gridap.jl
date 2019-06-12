module CellFieldsMocks

using Gridap
using Gridap.CellFieldsOperations: IterCellFieldLikeAndGradient
using Gridap.CellFieldsOperations: IndexCellFieldLikeAndGradient

using ..CellValuesMocks
using ..FieldsMocks

export IterCellFieldMock
export IndexCellFieldMock
export IterCellBasisMock
export IndexCellBasisMock

function IterCellFieldMock(D::Int,X::Type,l::Int)
  f = MockField(D,X)
  g = gradient(f)
  cf = TestIterCellValue(f,l)
  cg = TestIterCellValue(g,l)
  IterCellFieldLikeAndGradient(cf,cg)
end

function IndexCellFieldMock(D::Int,X::Type,l::Int)
  f = MockField(D,X)
  g = gradient(f)
  cf = TestIndexCellValue(f,l)
  cg = TestIndexCellValue(g,l)
  IndexCellFieldLikeAndGradient(cf,cg)
end

function IterCellBasisMock(D::Int,X::Type,l::Int)
  f = MockBasis(D,X)
  g = gradient(f)
  cf = TestIterCellValue(f,l)
  cg = TestIterCellValue(g,l)
  IterCellFieldLikeAndGradient(cf,cg)
end

function IndexCellBasisMock(D::Int,X::Type,l::Int)
  f = MockBasis(D,X)
  g = gradient(f)
  cf = TestIndexCellValue(f,l)
  cg = TestIndexCellValue(g,l)
  IndexCellFieldLikeAndGradient(cf,cg)
end

end # module
