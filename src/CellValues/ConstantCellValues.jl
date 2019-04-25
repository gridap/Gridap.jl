
# ConstantCellValue

struct ConstantCellValue{T} <: IndexCellValue{T,1}
  value::T
  length::Int
end

celldata(self::ConstantCellValue) = self.value

size(self::ConstantCellValue) = (self.length,)

IndexStyle(::Type{ConstantCellValue{T}} where T) = IndexLinear()

# ConstantCellArray

# struct ConstantCellArray{T,N} <: IndexCellArray{T,N,Array{T,N},1}
  # array::Array{T,N}
  # length::Int
# end

const ConstantCellArray{T,N} = ConstantCellValue{Array{T,N}}

const ConstantCellVector{T} = ConstantCellArray{T,1}

const ConstantCellMatrix{T} = ConstantCellArray{T,2}

cellsize(self::ConstantCellValue) = size(self.value)

function cellsum(self::ConstantCellArray{T,N};dim::Int) where {T,N}
  b = sum(self.value,dims=dim)
  s = cellsumsize(size(b),Val(dim))
  c = copy(reshape(b,s))
  ConstantCellValue(c,self.length)
end

function cellsum(self::ConstantCellArray{T,1};dim::Int) where T
  b = sum(self.value)
  ConstantCellValue(b,self.length)
end

function cellnewaxis(self::ConstantCellArray;dim::Int)
  s = [ v for v in size(self.value)]
  insert!(s,dim,1)
  shape = tuple(s...)
  c = copy(reshape(self.value,shape))
  ConstantCellValue(c,self.length)
end

# ConstantCellValue

getindex(self::ConstantCellValue,cell::Int) = celldata(self)

function (==)(a::ConstantCellValue,b::ConstantCellValue)
  celldata(a) != celldata(b) && return false
  length(a) != length(b) && return false
  return true
end

for op in (:+,:-,:*,:/,:(outer),:(inner))

  @eval begin
    function ($op)(a::ConstantCellValue,b::ConstantCellValue)
      @assert length(a) == length(b)
      c = broadcast($op,celldata(a),celldata(b))
      ConstantCellValue(c,length(a))
    end
  end

end

for op in (:+,:-,:(det),:(inv),:(meas))

  @eval begin
    function ($op)(a::ConstantCellValue)
      c = broadcast($op,celldata(a))
      ConstantCellValue(c,length(a))
    end
  end

end
