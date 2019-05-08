
function test_iter_cell_value(
  iter_cell_value::IterCellValue{T}, values::AbstractArray{T}) where T

  _test_iter_cell_value(iter_cell_value, values)
end

function test_index_cell_value(
  index_cell_value::IndexCellValue{T,N}, values::AbstractArray{T,N}) where {T,N}

  _test_index_cell_value(index_cell_value, values)
end

function test_iter_cell_array(
  iter_cell_array::IterCellArray{T,N},
  arrays::Array{Array{T,N}}) where {T,N}

  _test_iter_cell_value(iter_cell_array,arrays)

  cellsize(iter_cell_array) == maximum([ size(v) for v in arrays ])

end

function test_index_cell_array(
  index_cell_array::IndexCellArray{T,N},
  arrays::AbstractArray{Array{T,N}}) where {T,N}

  _test_index_cell_value(index_cell_array,arrays)

  cellsize(index_cell_array) == maximum([ size(v) for v in arrays ])

end

function _test_iter_cell_value(iter_cell_value, values)

  @test length(iter_cell_value) == length(values)

  i = 0
  for value in iter_cell_value
    i += 1
    @assert value == values[i]
  end

  @test i == length(values)

  @test eltype(iter_cell_value) == eltype(values)

  @test iter_cell_value == iter_cell_value

end

function _test_index_cell_value( index_cell_value, values)

  @test length(index_cell_value) == length(values)

  @test size(index_cell_value) == size(values)

  @test IndexStyle(index_cell_value) == IndexStyle(values)

  l = 0
  for i in eachindex(index_cell_value)
    @assert index_cell_value[i] == values[i]
    l +=1
  end

  @test l == length(values)

  @test index_cell_value == index_cell_value

end
