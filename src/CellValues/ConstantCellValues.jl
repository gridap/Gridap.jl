
# ConstantCellValue

struct ConstantCellValue{T} <: IndexCellValue{T,1}
  value::T
  length::Int
end

celldata(self::ConstantCellValue) = self.value

size(self::ConstantCellValue) = (self.length,)

IndexStyle(::Type{ConstantCellValue{T}} where T) = IndexLinear()

# ConstantCellArray

struct ConstantCellArray{T,N} <: IndexCellArray{T,N,Array{T,N},1}
  array::Array{T,N}
  length::Int
end

const ConstantCellVector{T} = ConstantCellArray{T,1}

const ConstantCellMatrix{T} = ConstantCellArray{T,2}

celldata(self::ConstantCellArray) = self.array

size(self::ConstantCellArray) = (self.length,)

IndexStyle(::Type{ConstantCellArray{T,N}} where {T,N}) = IndexLinear()

cellsize(self::ConstantCellArray) = size(self.array)

function cellsum(self::ConstantCellArray{T,N};dim::Int) where {T,N}
  b = sum(self.array,dims=dim)
  s = cellsumsize(size(b),Val(dim))
  c = copy(reshape(b,s))
  ConstantCellArray(c,self.length)
end

function cellsum(self::ConstantCellArray{T,1};dim::Int) where T
  b = sum(self.array)
  ConstantCellValue(b,self.length)
end

function cellnewaxis(self::ConstantCellArray;dim::Int)
  s = [ v for v in size(self.array)]
  insert!(s,dim,1)
  shape = tuple(s...)
  c = copy(reshape(self.array,shape))
  ConstantCellArray(c,self.length)
end

# ConstantCellData

const ConstantCellData{T} = Union{ConstantCellValue{T},ConstantCellArray{T}}

function ConstantCellData(a,l::Int)
  ConstantCellValue(a,l)
end

function ConstantCellData(a::AbstractArray,l::Int)
  ConstantCellArray(a,l)
end

celldata(::ConstantCellData) = @abstractmethod

getindex(self::ConstantCellData,cell::Int) = celldata(self)

function (==)(a::ConstantCellData,b::ConstantCellData)
  celldata(a) != celldata(b) && return false
  length(a) != length(b) && return false
  return true
end

for op in (:+,:-,:*,:/,:(outer),:(inner))

  @eval begin
    function ($op)(a::ConstantCellData,b::ConstantCellData)
      @assert length(a) == length(b)
      c = broadcast($op,celldata(a),celldata(b))
      ConstantCellData(c,length(a))
    end
  end

end

for op in (:+,:-,:(det),:(inv),:(meas))

  @eval begin
    function ($op)(a::ConstantCellData)
      c = broadcast($op,celldata(a))
      ConstantCellData(c,length(a))
    end
  end

end

