module CellValues

export CellValue
export CellArray
export CellMatrix
export CellVector

export IterCellValue
export IterCellArray
export IterCellMatrix
export IterCellVector

export IndexCellValue
export IndexCellArray
export IndexCellMatrix
export IndexCellVector

export apply
export cellsize
export celllength
export cellsum
export cellnewaxis
export cellmean

include("AbstractCellValues.jl")
using Numa.CellValues.AbstractCellValues

include("Operations.jl")
using Numa.CellValues.Operations

include("ConstantCellValues.jl")

include("Wrappers.jl")

include("Testers.jl")

end # module CellValues
