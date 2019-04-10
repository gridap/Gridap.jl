
struct FlexibleUnstructuredGrid{D,Z} <: Grid{D,Z}
  points::Vector{Point{D}}
  cells::Vector{Vector{Int}}
  ctypes::Vector{NTuple{Z,Int}}
end

points(self::FlexibleUnstructuredGrid) = CellValueFromArray(self.points)

cells(self::FlexibleUnstructuredGrid) = CellArrayFromArrayOfArrays(self.cells)

celltypes(self::FlexibleUnstructuredGrid) = CellValueFromArray(self.ctypes)

# Constructors

function FlexibleUnstructuredGrid(grid::CartesianGrid{D}) where D
  ps = Array{Point{D},1}(undef,(length(points(grid)),))
  for (i,xi) in enumerate(points(grid))
    ps[i] = xi
  end
  cs = [ Array{Int,1}(undef,(2^D,)) for i in 1:length(cells(grid)) ]
  for (i,ci) in enumerate(cells(grid))
    cs[i] .= ci
  end
  ts = [ ti for ti in celltypes(grid) ]
  FlexibleUnstructuredGrid(ps,cs,ts)
end

