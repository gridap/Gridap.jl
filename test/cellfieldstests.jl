
xe = Numa.DummyCellCoordinates2D(partition=(2,2))

basis = ShapeFunctionsScalarQua4()
cellbasis = Numa.CellBasisFromSingleInterpolation(basis)

phi = Numa.CellFieldFromInterpolation{2,Point{2},Float64,Point{2}}(cellbasis,xe)

quad = TensorProductQuadrature{2}(orders=[2,2])
cellquad = ConstantCellQuadrature(quad,length(xe))

q = coordinates(cellquad)

x = evaluate(phi,q)

for xcg in x
  for xg in xcg
  end
end

jac = gradient(phi)

jac_g = evaluate(jac,q)

using LinearAlgebra

function compute_domain_vol(jac_g)
  V = 0.0
  for jac_cg in jac_g
    for J in jac_cg
      V += det(J)
    end
  end
  V
end

V = compute_domain_vol(jac_g)
@test V ≈ 1.0


fun(x::Point{2}) = 2x[1]+x[2]

fun(::Type{Point{2}}) = Float64

g = compose(fun,phi)

@test isa(g,CellField{2,Float64})

gq = evaluate(g,q)

@test length(gq) == length(q)

@test maxsize(gq) == maxsize(q)

x = evaluate(phi,q)

for (fi,xi) in zip(gq,x)
  for (fii,xii) in zip(fi,xi)
    @assert fii == fun(xii)
  end
end


