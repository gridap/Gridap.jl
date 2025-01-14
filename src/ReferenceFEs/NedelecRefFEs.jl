
struct CurlConformity <: Conformity end

"""
    struct Nedelec <: PushforwardRefFE <: ReferenceFEName
"""
struct Nedelec <: PushforwardRefFE end

"""
    const nedelec = Nedelec()

Singleton of the [`Nedelec`](@ref) reference FE name.
"""
const nedelec = Nedelec()

Pushforward(::Type{<:Nedelec}) = CoVariantPiolaMap()

"""
    NedelecRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the curl of the  functions in this basis
is in the Q space of degree `order`.
"""
function NedelecRefFE(::Type{et},p::Polytope,order::Integer) where et
  D = num_dims(p)

  if is_n_cube(p)
    prebasis = QGradBasis(Monomial,Val(D),et,order) # Prebasis
    eb = MonomialBasis(Val(1),et,order)             # Edge basis
    fb = QGradBasis(Monomial,Val(D-1),et,order-1)   # Face basis
    cb = QCurlGradBasis(Monomial,Val(D),et,order-1) # Cell basis
  elseif is_simplex(p)
    prebasis = PGradBasis(Monomial,Val(D),et,order) # Prebasis
    eb = MonomialBasis(Val(1),et,order)             # Edge basis
    fb = MonomialBasis(Val(D-1),VectorValue{D-1,et},order-1,Polynomials._p_filter) # Face basis
    cb = MonomialBasis(Val(D),VectorValue{D,et},order-D+1,Polynomials._p_filter)   # Cell basis
  else
    @unreachable "Nedelec Reference FE only implemented for n-cubes and simplices"
  end

  function cmom(φ,μ,ds) # Cell moment function: σ_K(φ,μ) = ∫(φ⋅μ)dK
    Broadcasting(Operation(⋅))(φ,μ)
  end
  function fmom_HEX(φ,μ,ds) # Face moment function: σ_F(φ,μ) = ∫((φ×n)⋅μ)dF
    o = get_facet_orientations(ds.cpoly)[ds.face] # This is a hack to avoid a sign map
    n = o*get_facet_normal(ds)
    E = get_extension(ds)
    Eμ = Broadcasting(Operation(⋅))(E,μ) # We have to extend the basis to 3D
    φn = Broadcasting(Operation(×))(n,φ)
    Broadcasting(Operation(⋅))(φn,Eμ)
  end
  function fmom_TET(φ,μ,ds) # Face moment function: σ_F(φ,μ) = ∫((φ×n)⋅μ)dF
    E = get_extension(ds)
    Eμ = Broadcasting(Operation(⋅))(E,μ) # We have to extend the basis to 3D
    Broadcasting(Operation(⋅))(φ,Eμ)
  end
  function emom(φ,μ,ds) # Edge moment function: σ_E(φ,μ) = ∫((φ⋅t)*μ)dE
    t = get_edge_tangent(ds)
    φt = Broadcasting(Operation(⋅))(φ,t)
    Broadcasting(Operation(*))(φt,μ)
  end

  moments = Tuple[
    (get_dimrange(p,1),emom,eb), # Edge moments
  ]
  if (D == 3) && order > 0
    fmom = ifelse(is_n_cube(p),fmom_HEX,fmom_TET)
    push!(moments,(get_dimrange(p,D-1),fmom,fb)) # Face moments
  end
  if (is_n_cube(p) && order > 0) || (is_simplex(p) && order > D-2)
    push!(moments,(get_dimrange(p,D),cmom,cb))   # Cell moments
  end

  return MomentBasedReferenceFE(Nedelec(),p,prebasis,moments,CurlConformity())
end

function ReferenceFE(p::Polytope,::Nedelec, order)
  NedelecRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::Nedelec,::Type{T}, order) where T
  NedelecRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{Nedelec},sym::Symbol)
  hcurl = (:Hcurl,:HCurl)
  if sym == :L2
    L2Conformity()
  elseif sym in hcurl
    CurlConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a Nedelec reference FE.

    Possible values of conformity for this reference fe are $((:L2, hcurl...)).
    """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{Nedelec}, ::CurlConformity)
  reffe.face_dofs # For Nedelec, this member variable holds the face owned dofs
end

function get_face_dofs(reffe::GenericRefFE{Nedelec,Dc}) where Dc
  face_dofs = [Int[] for i in 1:num_faces(reffe)]
  face_own_dofs = get_face_own_dofs(reffe)
  p = get_polytope(reffe)
  for d = 1:Dc # Starting from edges, vertices do not own DoFs for Nedelec
    first_face = get_offset(p,d)
    nfaces     = num_faces(reffe,d)
    for face = first_face+1:first_face+nfaces
      for df = 1:d-1
        face_faces  = get_faces(p,d,df)
        first_cface = get_offset(p,df)
        for cface in face_faces[face-first_face]
          cface_own_dofs = face_own_dofs[first_cface+cface]
          for dof in cface_own_dofs
            push!(face_dofs[face],dof)
          end
        end
      end
      for dof in face_own_dofs[face]
        push!(face_dofs[face],dof)
      end
    end
  end
  face_dofs
end
