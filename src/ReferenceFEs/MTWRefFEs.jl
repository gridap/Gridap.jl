
struct MardalTaiWinther <: ReferenceFEName end

const mtw = MardalTaiWinther()

Pushforward(::Type{MardalTaiWinther}) = ContraVariantPiolaMap()

"""
    struct MardalTaiWinther <: ReferenceFEName end
    MardalTaiWintherRefFE(::Type{et},p::Polytope,order::Integer) where et

Mardal-Tai-Winther reference finite element.

References:

-2D: `A Robust Finite Element Method for Darcy-Stokes Flow`, Mardal, Tai and Winther (2002)
-3D: `A discrete de Rham complex with enhanced smoothness`, Tai and Winther (2006).
-mappings: `Transformations for Piola-mapped elements`, Aznaran, Farrell and Kirby (2022)
"""
function MardalTaiWintherRefFE(::Type{T},p::Polytope{D}) where {T,D}
  @check is_simplex(p) "MardalTaiWinther Reference FE only defined simplices."
  @check D in (2,3) "MardalTaiWinther Reference FE only defined in dimension 2 and 3."

  prebasis = _mtw_basis(Val(D),T)
  fb_n = FEEC_poly_basis(Val(D-1),T,1,0,:P⁻)
  # 2D: {1-x,x}       # in SEGMENT cartesian coords, so ≡ barycentric {λ1, λ2}
  # 3D: {1-x-y, x, y} # in TRI     cartesian coords, so ≡ barycentric {λ1, λ2, λ3}
  fb_t = FEEC_poly_basis(Val(D-1),T,1,1,:P⁻)
  # 2D: {1}
  # 3D: {(1-y, x), (y, 1-x), (-y, x)}

  function fmom_t(φ,μ,ds) # Face tangent moments: σ(φ,μ) = ∫((φ⋅μ)dK ≡ ∫((φ×n)⋅μ)dF
    E = get_extension(ds)
    Eμ = Broadcasting(Operation(⋅))(E,μ) # We have to extend the basis to 3D
    Broadcasting(Operation(⋅))(φ,Eμ)
  end
  function fmom_n(φ,μ,ds) # Face normal moments: σ(φ,μ) = ∫((φ·n)*μ)dF
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    Broadcasting(Operation(*))(φn,μ)
  end

  frange = get_dimrange(p,D-1)
  moments = [
    (frange,fmom_n, fb_n), # Face normal moments
    (frange,fmom_t, fb_t), # Face tangent moments
  ]

  return MomentBasedReferenceFE(mtw ,p,prebasis,moments,DivConformity())
end

function ReferenceFE(p::Polytope,::MardalTaiWinther,::Type{T}, order) where T
  @assert order == 3   "MardalTaiWinther Reference FE is only defined for order 3, got $order"
  MardalTaiWintherRefFE(T,p)
end


# pre-pre basis from DefElements (https://defelement.org/elements/examples/triangle-mardal-tai-winther-1.html)
function _mtw_basis(::Val{2},T)
  monoms = MonomialBasis(Val(2),VectorValue{2,T},3,Polynomials._p_filter)
  change = T[
#e₁ 1   x   x²   x³  y   xy   x²y y²   xy²  y³
#e₂   1   x   x²   x³  y    xy  x²y  y²   xy²  y³
    1 0 0 0 0 0  0 0 0 0 0  0 0 0 0  0 0  0 0  0; # (1, 0)
    0 0 1 0 0 0  0 0 0 0 0  0 0 0 0  0 0  0 0  0; # (x, 0)
    0 0 0 0 0 0  0 0 1 0 0  0 0 0 0  0 0  0 0  0; # (y, 0)
    0 1 0 0 0 0  0 0 0 0 0  0 0 0 0  0 0  0 0  0; # (0, 1)
    0 0 0 1 0 0  0 0 0 0 0  0 0 0 0  0 0  0 0  0; # (0, x)
    0 0 0 0 0 0  0 0 0 1 0  0 0 0 0  0 0  0 0  0; # (0, y)
    0 0 0 0 1 0  0 0 0 0 2 -2 0 0 0 -1 0  0 0  0; # (x² + 2xy, -2xy -y²)
    0 0 0 0 2 0 -1 0 0 0 0 -4 0 3 0  0 3  0 0 -1; # (-x³ + 2x² + 3xy², 3x²y - 4xy -y³)
    0 0 0 0 1 0  0 0 0 0 0 -2 2 0 0  0 3 -2 0 -1] # (2x²y + x² + 3xy², -2xy² - 2xy -y³)
#∅            X    X              X         X   ] could optimise by filtering these terms from `monoms` using a CompWiseTensorPolyBasis.
  linear_combination(permutedims(change), monoms)
end

# pre-pre basis from DefElements (https://defelement.org/elements/examples/tetrahedron-mardal-tai-winther-1.html)
function _mtw_basis(::Val{3},::Type{T}) where T
  poly_to_comp_to_coefsAndExponents = Vector{Vector{Tuple{T,Tuple{Int,Int,Int}}}}[
    # (1, 0, 0) (0, 1, 0) (0, 0, 1)
    [ [(1, (0, 0, 0))], [], [] ],
    [ [], [(1, (0, 0, 0)) ], []],
    [ [], [], [(1, (0, 0, 0))] ],
    # (x, 0, 0) (0, x, 0) (0, 0, x)
    [ [(1, (1, 0, 0))], [], [] ],
    [ [], [(1, (1, 0, 0)) ], []],
    [ [], [], [(1, (1, 0, 0))] ],
    # (y, 0, 0) (0, y, 0) (0, 0, y)
    [ [(1, (0, 1, 0))], [], [] ],
    [ [], [(1, (0, 1, 0)) ], []],
    [ [], [], [(1, (0, 1, 0))] ],
    # (z, 0, 0) (0, z, 0) (0, 0, z)
    [ [(1, (0, 0, 1))], [], [] ],
    [ [], [(1, (0, 0, 1)) ], []],
    [ [], [], [(1, (0, 0, 1))] ],
    # (0, -x²y -xy² -2xyz +xy,  x²z +2xyz +xz² -xz)
    [ [], [(-1, (2, 1, 0)), (-1, (1, 2, 0)), (-2, (1, 1, 1)), (1, (1, 1, 0))], [(1, (2, 0, 1)), (2, (1, 1, 1)), (1, (1, 0, 2)), (-1, (1, 0, 1)) ] ],
    # (x²y +xy² +2xyz -xy,  0,  -y²z -2xyz -yz² +yz)
    [ [(1, (2, 1, 0)), (1, (1, 2, 0)), (2, (1, 1, 1)), (-1, (1, 1, 0))], [], [(-1, (0, 2, 1)), (-2, (1, 1, 1)), (-1, (0, 1, 2)), (1, (0, 1, 1)) ] ],
    # (-x²z -2xyz -xz² +xz,  y²z +2xyz +yz² -yz,  0)
    [ [(-1, (2, 0, 1)), (-2, (1, 1, 1)), (-1, (1, 0, 2)), (1, (1, 0, 1)) ], [(1, (0, 2, 1)), (2, (1, 1, 1)), (1, (0, 1, 2)), (-1, (0, 1, 1)) ], [] ],
    # (0,  -x³y -x²y² -2x²yz +x²y,  x³z +2x²yz +x²z² -x²z)
    [ [], [(-1, (3, 1, 0)), (-1, (2, 2, 0)), (-2, (2, 1, 1)), (1, (2, 1, 0))], [(1, (3, 0, 1)), (2, (2, 1, 1)), (1, (2, 0, 2)), (-1, (2, 0, 1)) ] ],
    # (x³y +x²y² +2x²yz -x²y,  0,  -3x²yz -2xy²z -2xyz² +2xyz)
    [ [(1, (3, 1, 0)), (1, (2, 2, 0)), (2, (2, 1, 1)), (-1, (2, 1, 0))], [], [(-3, (2, 1, 1)), (-2, (1, 2, 1)), (-2, (1, 1, 2)), (2, (1, 1, 1)) ] ],
    # (-x³z -2x²yz -x²z² +x²z,  3x²yz +2xy²z +2xyz² -2xyz,  0)
    [ [(-1, (3, 0, 1)), (-2, (2, 1, 1)), (-1, (2, 0, 2)), (1, (2, 0, 1)) ], [(3, (2, 1, 1)), (2, (1, 2, 1)), (2, (1, 1, 2)), (-2, (1, 1, 1)) ], [] ],
    # (0,  -x²y² -xy³ -2xy²z +xy²,  2x²yz +3xy²z +2xyz² -2xyz)
    [ [], [(-1, (2, 2, 0)), (-1, (1, 3, 0)), (-2, (1, 2, 1)), (1, (1, 2, 0))], [(2, (2, 1, 1)), (3, (1, 2, 1)), (2, (1, 1, 2)), (-2, (1, 1, 1)) ] ],
    # (x²y² +xy³ +2xy²z -xy²,  0,  -2xy²z -y³z -y²z² +y²z)
    [ [(1, (2, 2, 0)), (1, (1, 3, 0)), (2, (1, 2, 1)), (-1, (1, 2, 0))], [], [(-2, (1, 2, 1)), (-1, (0, 3, 1)), (-1, (0, 2, 2)), (1, (0, 2, 1)) ] ],
    # (-2x²yz -3xy²z -2xyz² +2xyz, 2xy²z +y³z +y²z² -y²z,  0)
    [ [(-2, (2, 1, 1)), (-3, (1, 2, 1)), (-2, (1, 1, 2)), (2, (1, 1, 1)) ], [(2, (1, 2, 1)), (1, (0, 3, 1)), (1, (0, 2, 2)), (-1, (0, 2, 1)) ], [] ],
    # (0,  -2x²yz -2xy²z -3xyz² +2xyz,  x²z² +xz³ +2xyz² -xz²)
    [ [], [(-2, (2, 1, 1)), (-2, (1, 2, 1)), (-3, (1, 1, 2)), (2, (1, 1, 1)) ], [(1, (2, 0, 2)), (1, (1, 0, 3)), (2, (1, 1, 2)), (-1, (1, 0, 2)) ] ],
    # (2x²yz +2xy²z +3xyz² -2xyz,  0,  -y²z² -yz³ -2xyz² +yz²)
    [ [(2, (2, 1, 1)), (2, (1, 2, 1)), (3, (1, 1, 2)), (-2, (1, 1, 1)) ], [], [(-1, (0, 2, 2)), (-1, (0, 1, 3)), (-2, (1, 1, 2)), (1, (0, 1, 2)) ] ],
    # (-x²z² -xz³ -2xyz² +xz²,  y²z² +yz³ +2xyz² -yz²,  0)
    [ [(-1, (2, 0, 2)), (-1, (1, 0, 3)), (-2, (1, 1, 2)), (1, (1, 0, 2)) ], [(1, (0, 2, 2)), (1, (0, 1, 3)), (2, (1, 1, 2)), (-1, (0, 1, 2)) ], [] ],
  ]

  comp_terms = [ CartesianIndex{3}[] for _ in 1:3 ]
  for comp_to_coefsAndExponents in poly_to_comp_to_coefsAndExponents
    for (l, coefsAndExponents) in enumerate(comp_to_coefsAndExponents)
      for (_, exponents) in coefsAndExponents
        term = CartesianIndex(exponents .+ 1)
        term in comp_terms[l] || push!(comp_terms[l], CartesianIndex(term))
      end
    end
  end

  orders = [3 3 3; 3 3 3; 3 3 3;]
  poly_prebasis = CompWiseTensorPolyBasis{3}(Monomial, VectorValue{3,T}, orders, comp_terms)

  num_polys = length(poly_to_comp_to_coefsAndExponents)
  num_pb_polys = mapreduce(length, +, comp_terms)
  change = zeros(T, (num_pb_polys, num_polys))
  comps_num_term = length.(comp_terms)

  # need to loop again because we don't know the final number of terms per comp
  # initially, so the size of change was unknown
  for (poly, comp_to_coefsAndExponents) in enumerate(poly_to_comp_to_coefsAndExponents)
    for (l,coefsAndExponents) in enumerate(comp_to_coefsAndExponents)
      offset = sum(comps_num_term[1:l-1])
      for (coef, exponents) in coefsAndExponents
        term = CartesianIndex(exponents .+ 1)
        i = findfirst(==(term), comp_terms[l])
        change[offset+i,poly] = coef
      end
    end
  end

  linear_combination(change, poly_prebasis)
end

