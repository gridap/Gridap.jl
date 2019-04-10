
struct UnstructuredGrid{D,Z} <: Grid{D,Z}
  points::Vector{Point{D}}
  cells_data::Vector{Int}
  cells_ptrs::Vector{Int}
  ctypes::Vector{NTuple{Z,Int}}
end

points(self::UnstructuredGrid) = CellValueFromArray(self.points)

cells(self::UnstructuredGrid) = CellVectorFromDataAndPtrs(self.cells_data,self.cells_ptrs)

celltypes(self::UnstructuredGrid) = CellValueFromArray(self.ctypes)

struct FlexibleUnstructuredGrid{D,Z} <: Grid{D,Z}
  points::Vector{Point{D}}
  cells::Vector{Vector{Int}}
  ctypes::Vector{NTuple{Z,Int}}
end

points(self::FlexibleUnstructuredGrid) = CellValueFromArray(self.points)

cells(self::FlexibleUnstructuredGrid) = CellArrayFromArrayOfArrays(self.cells)

celltypes(self::FlexibleUnstructuredGrid) = CellValueFromArray(self.ctypes)

# Constructors from CartesianGrid (mainly for testing purposes)

function UnstructuredGrid(grid::CartesianGrid{D}) where D
  ps = compute_points(grid)
  ts = compute_celltypes(grid)
  data, ptrs = compute_cells(grid,UnstructuredGrid{D,D})
  UnstructuredGrid(ps,data,ptrs,ts)
end

function FlexibleUnstructuredGrid(grid::CartesianGrid{D}) where D
  ps = compute_points(grid)
  ts = compute_celltypes(grid)
  cs = compute_cells(grid,FlexibleUnstructuredGrid{D,D})
  FlexibleUnstructuredGrid(ps,cs,ts)
end

function compute_points(grid::CartesianGrid{D}) where D
  ps = Array{Point{D},1}(undef,(length(points(grid)),))
  for (i,xi) in enumerate(points(grid))
    ps[i] = xi
  end
  ps
end

function compute_celltypes(grid::CartesianGrid)
  [ ti for ti in celltypes(grid) ]
end

function compute_cells(grid::CartesianGrid{D},::Type{FlexibleUnstructuredGrid{D,D}}) where D
  cs = [ Array{Int,1}(undef,(2^D,)) for i in 1:length(cells(grid)) ]
  for (i,ci) in enumerate(cells(grid))
    cs[i] .= ci
  end
  cs
end

function compute_cells(grid::CartesianGrid{D},::Type{UnstructuredGrid{D,D}}) where D
  ptrs = fill(2^D,(length(cells(grid))+1,))
  length_to_ptrs!(ptrs)
  data = zeros(Int,ptrs[end]-1)
  fill_cell_data!(data,cells(grid))
  (data, ptrs)
end

function fill_cell_data!(data,cells)
  k = 1
  for v in cells
    for vi in v
      @inbounds data[k] = vi
      k +=1
    end
  end
end

