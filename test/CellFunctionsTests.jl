module CellFieldsTests

using Test
using Numa.FieldValues
using Numa.CellArrays
using Numa.CellFunctions
using Numa.Quadratures
using Numa.CellQuadratures

l = 10

include("CellFunctionsTestsMocks.jl")

@testset "InnerFieldValuesFieldValues" begin

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

@testset "InnerBasisValuesFieldValues" begin

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

@testset "InnerBasisValuesBasisValues" begin

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

@testset "ExpandBasisValuesFieldValues" begin

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

@testset "CellBasisFromSingleInterpolation" begin

  include("PolynomialsTestsMocks.jl")

  using Numa.CellFunctions: CellBasisValuesFromSingleInterpolation

  l = 10

  refquad = TensorProductQuadrature(orders=(5,4))
  refpoints = coordinates(refquad)

  quad = ConstantCellQuadrature(refquad,l)
  points = coordinates(quad)

  refbasis = ShapeFunctionsScalarQua4()

  refvals = evaluate(refbasis,refpoints)

  vals = CellBasisValuesFromSingleInterpolation(refbasis,points)

  @test isa(vals,CellBasisValues{Float64})

  for refvals2 in vals
    @assert refvals2 == refvals
  end

  basis = CellBasisFromSingleInterpolation(refbasis)

  @test isa(basis,CellBasis{2,Float64})

  vals = evaluate(basis,points)

  @test isa(vals,CellBasisValues{Float64})

  @test isa(vals,ConstantCellArray{Float64,2})

  for refvals2 in vals
    @assert refvals2 == refvals
  end

  vals = evaluate(basis,vfv)

  @test isa(vals,CellBasisValues{Float64})

  @test isa(vals,CellBasisValuesFromSingleInterpolation)

  refbasisgrad = gradient(refbasis)

  refvalsgrad = evaluate(refbasisgrad,refpoints)

  basisgrad = gradient(basis)

  valsgrad = evaluate(basisgrad,points)

  @test isa(valsgrad,CellBasisValues{VectorValue{2}})

  @test isa(valsgrad,ConstantCellArray{VectorValue{2},2})

  for refvalsgrad2 in valsgrad
    @assert refvalsgrad2 == refvalsgrad
  end

  valsgrad = CellBasisValuesFromSingleInterpolation(refbasisgrad,points)

  @test isa(valsgrad,CellBasisValues{VectorValue{2}})

  for refvalsgrad2 in valsgrad
    @assert refvalsgrad2 == refvalsgrad
  end

end

@testset "CellBasisOps" begin

  sb = ScalarBasisMock()

  vals = evaluate(sb + sb,vfv)

  @test isa(sb + sb,CellBasis{2,Float64})

  for a in vals
    @assert a == sbva + sbva
  end

  vals = evaluate(sb - sb,vfv)

  @test isa(sb - sb,CellBasis{2,Float64})

  for a in vals
    @assert a == sbva - sbva
  end

  vals = evaluate(sb * sb,vfv)

  @test isa(sb * sb,CellBasis{2,Float64})

  for a in vals
    @assert a == sbva .* sbva
  end

  vals = evaluate(sb / sb,vfv)

  @test isa(sb / sb,CellBasis{2,Float64})

  for a in vals
    @assert a == sbva ./ sbva
  end

  vals = evaluate(inner(sb,sb),vfv)

  @test isa(vals,CellArray{Float64,3})

  sf = ScalarFieldMock()

  vals = evaluate(inner(sf,sf),vfv)

  @test isa(vals,CellArray{Float64,1})

  vals = evaluate(inner(sb,sf),vfv)

  @test isa(vals,CellArray{Float64,2})

end

@testset "CellBasisWithGeomap" begin

  include("IntegrationMeshesTestsMocks.jl")

  imesh = DummyIntegrationMesh2D(partition=(3,3))
  refquad = TensorProductQuadrature(orders=(0,0))
  meshcoords = cellcoordinates(imesh)
  quad = ConstantCellQuadrature(refquad,length(meshcoords))
  points = coordinates(quad)
  phi = geomap(imesh)

  basis = cellbasis(imesh)

  physbasis = attachgeomap(basis,phi)

  vals = evaluate(physbasis,points)

  @test isa(vals,ConstantCellArray{Float64,2})

  physbasisgrad = gradient(physbasis)

  valsgrad = evaluate(physbasisgrad,points)

  tv1 = VectorValue(-1.5, -1.5)
  tv2 = VectorValue(1.5, -1.5)
  tv3 = VectorValue(-1.5, 1.5)
  tv4 = VectorValue(1.5, 1.5)
  
  valsgradref = reshape([tv1, tv2, tv3, tv4],(4,1))

  for v in valsgrad
    @assert v ≈ valsgradref
  end

end

end  # module CellFieldsTests
