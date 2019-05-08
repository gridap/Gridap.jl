struct IndexCellValueByGlobalAppend{T,V<:IndexCellValue{T,1},W<:IndexCellValue{T,1}} <: IndexCellValue{T,1}
  cvs1::V
  cvs2::W
  offset::Int
end

@propagate_inbounds function getindex(self::IndexCellValueByGlobalAppend,cell::Int)
  if cell <= self.offset
    return self.cvs1[cell]
  else
    return self.cvs2[cell-self.offset]
  end
end

size(self::IndexCellValueByGlobalAppend) = (length(self.cvs1) + length(self.cvs2),)

IndexStyle(::Type{IndexCellValueByGlobalAppend{T,N,V}}) where {T,N,V} = IndexStyle(V)

function IndexCellValueByGlobalAppend(cvs1::IndexCellValue{T,1}, cvs2::IndexCellValue{T,1}) where {T}
  offset = length(cvs1)
  IndexCellValueByGlobalAppend(cvs1, cvs2, offset)
end

function IndexCellValueByGlobalAppend(cvs::Vararg{IndexCellValue{T,1},N}) where {T,N}
  aux = cvs[1]
  for i in 2:N
    aux = IndexCellValueByGlobalAppend(aux, cvs[i])
  end
  return aux
end
