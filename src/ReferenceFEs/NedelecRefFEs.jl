"""
    struct Nedelec{kind} <: ReferenceFEName
"""
struct Nedelec{kind} <: ReferenceFEName
  Nedelec{1}() = new{1}()
  Nedelec{2}() = new{2}()
end

"""
    const nedelec = Nedelec{1}()
    const nedelec1 = nedelec

Singleton of the first kind of [`Nedelec`](@ref) reference FE name.
"""
const nedelec = Nedelec{1}()
const nedelec1 = nedelec

"""
    const nedelec2 = Nedelec{2}()

Singleton of the second kind of [`Nedelec`](@ref) reference FE name.
"""
const nedelec2 = Nedelec{2}()

Pushforward(::Type{<:Nedelec}) = CoVariantPiolaMap()

"""
    NedelecRefFE(::Type{T}, p::Polytope, order::Integer; sh_is_pb=true)

The `order` argument has the following meaning: the curl of the  functions in
this basis is in the ℙ/ℚ space of degree `order`. `T` is the type of scalar components.

`sh_is_pb` is only used if `p` is a simplex.
"""
function NedelecRefFE(::Type{T},p::Polytope,order::Integer; sh_is_pb=true, kind::Int=1) where T
  D = num_dims(p)
  rotate_90 = D==2

  if is_n_cube(p)
    PT = Legendre # Could be a wargs, any basis works
    @check kind == 1 "Nedelec reference elements of the second kind are only defined on simplices"
    prebasis =     FEEC_poly_basis(Val(D),T,order+1,1,:Q⁻,PT) # Q⁻ᵣΛ¹(□ᴰ), r = order+1
    eb =           FEEC_poly_basis(Val(1),T,order,0,  :Q⁻,PT)                      # Edge basis  Q⁻ᵨΛ⁰(□¹),  ρ = r-1
    fb = order>0 ? FEEC_poly_basis(Val(2),T,order,1,  :Q⁻,PT)            : nothing # Facet basis Q⁻ᵨΛ¹(□²),  ρ = r-1 (only D=3)
    cb = order>0 ? FEEC_poly_basis(Val(D),T,order,D-1,:Q⁻,PT; rotate_90) : nothing # Cell basis  Q⁻ᵨΛᴰ⁻¹(□ᴰ),ρ = r-1
  elseif is_simplex(p)
    PT = Bernstein # Could be a Kwargs, any basis works
    if kind == 1
      prebasis =       FEEC_poly_basis(Val(D),T,order+1,1,    :P⁻,PT) # P⁻ᵣΛ¹(△ᴰ), r = order+1
      eb =             FEEC_poly_basis(Val(1),T,order,0,      :P ,PT)                           # Edge basis  PᵨΛ⁰(△¹),  ρ = r-1
      fb = order>0 ?   FEEC_poly_basis(Val(2),T,order-1,1,    :P ,PT; rotate_90=true) : nothing # Facet basis PᵨΛ¹(△²),  ρ = r-2 (only D=3)
      cb = order>D-2 ? FEEC_poly_basis(Val(D),T,order-D+1,D-1,:P ,PT; rotate_90)      : nothing # Cell basis  PᵨΛᴰ⁻¹(△ᴰ),ρ = r-D
    else
      @check order > 0 "the lowest order of Nedelec elements of second kind is 1"
      prebasis =     FEEC_poly_basis(Val(D),T,order,1,      :P ,PT) # PᵣΛ¹(△ᴰ), r = order
      eb =           FEEC_poly_basis(Val(1),T,order,0,      :P⁻,PT)                           # Edge basis  P⁻ᵨΛ⁰(△¹),  ρ = r
      fb = order>1 ? FEEC_poly_basis(Val(2),T,order-1,1,    :P⁻,PT; rotate_90=true) : nothing # Facet basis P⁻ᵨΛ¹(△²),  ρ = r-1 (only D=3)
      cb = order≥D ? FEEC_poly_basis(Val(D),T,order-D+1,D-1,:P⁻,PT; rotate_90)      : nothing # Cell basis  P⁻ᵨΛᴰ⁻¹(△ᴰ),ρ = r-D+1
    end
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
  if (D == 3) && order > kind-1
    fmom = ifelse(is_n_cube(p),fmom_HEX,fmom_TET)
    push!(moments,(get_dimrange(p,D-1),fmom,fb)) # Face moments
  end
  if (is_n_cube(p) && order > 0) || (is_simplex(p) && order > D+kind-3)
    push!(moments,(get_dimrange(p,D),cmom,cb))   # Cell moments
  end

  sh_is_pb = sh_is_pb && is_simplex(p)
  return MomentBasedReferenceFE(Nedelec{kind}(),p,prebasis,moments,CurlConformity(); sh_is_pb)
end

function ReferenceFE(p::Polytope,::Nedelec{K},::Type{T}, order; kwargs...) where {K, T}
  NedelecRefFE(T,p,order; kind=K, kwargs...)
end

function Conformity(reffe::GenericRefFE{<:Nedelec},sym::Symbol)
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

function get_face_own_dofs(reffe::GenericRefFE{<:Nedelec}, ::CurlConformity)
  reffe.face_dofs # For Nedelec, this member variable holds the face owned dofs
end

function get_face_dofs(reffe::GenericRefFE{<:Nedelec,Dc}) where Dc
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
