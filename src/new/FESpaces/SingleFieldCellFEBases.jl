
abstract type SingleFieldCellFEBasis <: GridapType end

"""
"""
function get_array(cb::SingleFieldCellFEBasis)
  @abstractmethod
end

"""
"""
function get_cell_map(cb::SingleFieldCellFEBasis)
  @abstractmethod
end

struct CellShapeFunsWithMap{A,B} <: SingleFieldCellFEBasis
  cell_shapefuns::A
  cell_map::B
  @doc """
  """
  function CellShapeFunsWithMap(
    cell_shapefuns::AbstractArray{<:Field},
    cell_map::AbstractArray{<:Field})

    A = typeof(cell_shapefuns)
    B = typeof(cell_map)
    new{A,B}(cell_shapefuns,cell_map)
  end
end

function get_array(cb::CellShapeFunsWithMap)
  cb.cell_shapefuns
end

function get_cell_map(cb::CellShapeFunsWithMap)
  cb.cell_map
end

