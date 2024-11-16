struct BDM <: DivConforming end

const bdm = BDM()

"""
BDMRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the P space of degree `order-1`.

"""
function BDMRefFE(::Type{et},p::Polytope,order::Integer) where et
  D = num_dims(p)
  vet = VectorValue{D,et}

  if is_simplex(p)
    prebasis = MonomialBasis(vet,p,order)
    fb = MonomialBasis{D-1}(et,order,Polynomials._p_filter)
    cb = Polynomials.NedelecPrebasisOnSimplex{D}(order-2)
  else
    @notimplemented "BDM Reference FE only available for simplices"
  end

  function cmom(φ,μ,ds) # Cell moment function: σ_K(φ,μ) = ∫(φ·μ)dK
    Broadcasting(Operation(⋅))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function : σ_F(φ,μ) = ∫((φ·n)*μ)dF
    n = get_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    Broadcasting(Operation(*))(φn,μ)
  end
  moments = [
    (get_dimrange(p,D-1),fmom,fb), # Face moments
    (get_dimrange(p,D),cmom,cb)    # Cell moments
  ]

  return MomentBasedReferenceFE(BDM(),p,prebasis,moments,DivConformity())
end

function ReferenceFE(p::Polytope,::BDM, order)
  BDMRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::BDM,::Type{T}, order) where T
  BDMRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{BDM},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a BDM reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
      """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{BDM}, conf::DivConformity)
  get_face_dofs(reffe)
end
