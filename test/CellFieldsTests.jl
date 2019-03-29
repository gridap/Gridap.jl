module CellFieldsTests

using Test
using Numa.FieldValues
using Numa.CellArrays
using Numa.CellFields

import Numa.CellArrays: cellsize

struct CellArrayMockup{T,N} <: IndexableCellArray{T,N}
  a::Array{T,N}
  l::Int
end

function Base.getindex(self::CellArrayMockup,cell::Int)
  @assert 1 <= cell
  @assert cell <= length(self)
  self.a
end

Base.length(self::CellArrayMockup) = self.l

cellsize(self::CellArrayMockup) = size(self.a)

const CellFieldValuesMockup{T} = CellArrayMockup{T,1} where T <:FieldValue
const CellBasisValuesMockup{T} = CellArrayMockup{T,2} where T <:FieldValue

l = 10

sfva = [1.0,2.3,3.1,3.2]
sbva = [sfva'; sfva'; sfva']
sfv = CellFieldValuesMockup{Float64}(sfva,l)
sbv = CellBasisValuesMockup{Float64}(sbva,l)

vv = VectorValue(2.0,2.1)
vfva = [vv,2*vv,-1*vv,0.3*vv]
vbva = Array{VectorValue{2},2}(undef,(3,4))
for j in 1:4
  for i in 1:3
    vbva[i,j] = (i*j)*vv
  end
end
vfv = CellFieldValuesMockup{VectorValue{2}}(vfva,l)
vbv = CellBasisValuesMockup{VectorValue{2}}(vbva,l)

tv = TensorValue(2.0,2.1,3.4,4.1)
tfva = [tv,2*tv,-1*tv,0.3*tv]
tbva = Array{TensorValue{2,4},2}(undef,(3,4))
for j in 1:4
  for i in 1:3
    tbva[i,j] = (i*j)*tv
  end
end
tfv = CellFieldValuesMockup{TensorValue{2,4}}(tfva,l)
tbv = CellBasisValuesMockup{TensorValue{2,4}}(tbva,l)

@testset "InnerFieldField" begin

  siff = inner(sfv,sfv)
  siffa = [ inner(sfva[i],sfva[i]) for i in 1:4 ]

  @test length(siff) == l
  @test cellsize(siff) == size(sfva)
  for a in siff
    @assert a == siffa
  end

  viff = inner(vfv,vfv)
  viffa = [ inner(vfva[i],vfva[i]) for i in 1:4 ]

  @test length(viff) == l
  @test cellsize(viff) == size(vfva)
  for a in viff
    @assert a == viffa
  end

  tiff = inner(tfv,tfv)
  tiffa = [ inner(tfva[i],tfva[i]) for i in 1:4 ]

  @test length(tiff) == l
  @test cellsize(tiff) == size(tfva)
  for a in tiff
    @assert a == tiffa
  end

end

@testset "InnerBasisField" begin

  siff = inner(sbv,sfv)
  siffa = zeros(3,4)
  for j in 1:4
    for i in 1:3
      siffa[i,j] = inner(sbva[i,j],sfva[j])
    end
  end

  @test length(siff) == l
  @test cellsize(siff) == size(siffa)
  for a in siff
    @assert a == siffa
  end

  viff = inner(vbv,vfv)
  viffa = zeros(3,4)
  for j in 1:4
    for i in 1:3
      viffa[i,j] = inner(vbva[i,j],vfva[j])
    end
  end

  @test length(viff) == l
  @test cellsize(viff) == size(viffa)
  for a in viff
    @assert a == viffa
  end

  tiff = inner(tbv,tfv)
  tiffa = zeros(3,4)
  for j in 1:4
    for i in 1:3
      tiffa[i,j] = inner(tbva[i,j],tfva[j])
    end
  end

  @test length(tiff) == l
  @test cellsize(tiff) == size(tiffa)
  for a in tiff
    @assert a == tiffa
  end

end

@testset "InnerBasisBasis" begin

  siff = inner(sbv,sbv)
  siffa = zeros(3,3,4)
  for j in 1:4
    for i in 1:3
      for k in 1:3
        siffa[k,i,j] = inner(sbva[k,j],sbva[i,j])
      end
    end
  end

  @test length(siff) == l
  @test cellsize(siff) == size(siffa)
  for a in siff
    @assert a == siffa
  end

  viff = inner(vbv,vbv)
  viffa = zeros(3,3,4)
  for j in 1:4
    for i in 1:3
      for k in 1:3
        viffa[k,i,j] = inner(vbva[k,j],vbva[i,j])
      end
    end
  end

  @test length(viff) == l
  @test cellsize(viff) == size(viffa)
  for a in viff
    @assert a == viffa
  end

  tiff = inner(tbv,tbv)
  tiffa = zeros(3,3,4)
  for j in 1:4
    for i in 1:3
      for k in 1:3
        tiffa[k,i,j] = inner(tbva[k,j],tbva[i,j])
      end
    end
  end

  @test length(tiff) == l
  @test cellsize(tiff) == size(tiffa)
  for a in tiff
    @assert a == tiffa
  end

end

@testset "Expand" begin

  sfv2 = CellFieldValuesMockup{Float64}(sfva[1:3],l)

  sexpand = expand(sbv,sfv2)
  sexpanda = Array{Float64,1}(undef,(4,))
  for j in 1:4
    sexpanda[j] = 0.0
    for i in 1:3
      sexpanda[j] += sbva[i,j]*sfva[i]
    end
  end

  @test length(sexpand) == l
  @test cellsize(sexpand) == size(sexpanda)
  for a in sexpand
    @assert a == sexpanda
  end

  vfv2 = CellFieldValuesMockup{VectorValue{2}}(vfva[1:3],l)

  vexpand = expand(sbv,vfv2)
  vexpanda = Array{VectorValue{2},1}(undef,(4,))
  for j in 1:4
    vexpanda[j] = zero(VectorValue{2})
    for i in 1:3
      vexpanda[j] += sbva[i,j]*vfva[i]
    end
  end

  @test length(vexpand) == l
  @test cellsize(vexpand) == size(vexpanda)
  for a in vexpand
    @assert a == vexpanda
  end

  vexpand = expand(vbv,sfv2)
  vexpanda = Array{VectorValue{2},1}(undef,(4,))
  for j in 1:4
    vexpanda[j] = zero(VectorValue{2})
    for i in 1:3
      vexpanda[j] += vbva[i,j]*sfva[i]
    end
  end

  @test length(vexpand) == l
  @test cellsize(vexpand) == size(vexpanda)
  for a in vexpand
    @assert a == vexpanda
  end

end

end  # module CellFieldsTests
