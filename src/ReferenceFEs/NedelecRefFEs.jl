"""
    struct Nedelec <: ReferenceFEName
"""
struct Nedelec <: ReferenceFEName end

"""
    const nedelec = Nedelec()

Singleton of the [`Nedelec`](@ref) reference FE name.
"""
const nedelec = Nedelec()

Pushforward(::Type{Nedelec}) = CoVariantPiolaMap()

"""
    NedelecRefFE(::Type{T}, p::Polytope, order::Integer)

The `order` argument has the following meaning: the curl of the  functions in
this basis is in the ℙ/ℚ space of degree `order`. `T` is the type of scalar components.
"""
function NedelecRefFE(::Type{T},p::Polytope,order::Integer) where T
  D = num_dims(p)
  rotate_90 = D==2

  if is_n_cube(p)
    PT = Legendre # Could be a Kwargs, any hierarchical basis works
    @check PT ≠ Bernstein # broken for Bernstein in prebasis or cb, might be an issue of ordering of the basis polynomials
    prebasis =     FEEC_poly_basis(Val(D),T,order+1,1,:Q⁻,PT;) # Q⁻ᵣΛ¹(□ᴰ), r = order+1
    eb =           FEEC_poly_basis(Val(1),T,order,0,  :Q⁻,PT;)                     # Edge basis  Q⁻ᵨΛ⁰(□¹), ρ = r-1
    fb = order>0 ? FEEC_poly_basis(Val(2),T,order,1,  :Q⁻,PT; rotate_90) : nothing # Facet basis Q⁻ᵨΛ¹(□²), ρ = r-1 (only D=3)
    cb = order>0 ? FEEC_poly_basis(Val(D),T,order,D-1,:Q⁻,PT; rotate_90) : nothing # Cell basis  Q⁻ᵨΛ¹(□ᴰ), ρ = r-1
  elseif is_simplex(p)
    PT = Bernstein # Could be a Kwargs, any basis works
    prebasis =       FEEC_poly_basis(Val(D),T,order+1,1,  :P⁻,PT) # P⁻ᵣΛ¹(△ᴰ), r = order+1
    eb =             FEEC_poly_basis(Val(1),T,order,0,    :P ,PT)                      # Edge basis  PᵨΛ⁰(△¹), ρ = r-1
    fb = order>0 ?   FEEC_poly_basis(Val(2),T,order-1,1,  :P ,PT; rotate_90) : nothing # Facet basis PᵨΛ¹(△²), ρ = r-2 (only D=3)
    cb = order>D-2 ? FEEC_poly_basis(Val(D),T,order-D+1,1,:P ,PT; rotate_90) : nothing # Cell basis  PᵨΛ¹(△ᴰ), ρ = r-D
  else
    @notimplemented "Nedelec Reference FE only implemented for n-cubes and simplices"
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
