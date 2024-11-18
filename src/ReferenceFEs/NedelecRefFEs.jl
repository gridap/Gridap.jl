
struct CurlConformity <: Conformity end

struct Nedelec <: ReferenceFEName end
const nedelec = Nedelec()

"""
    NedelecRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the curl of the  functions in this basis
is in the Q space of degree `order`.
"""
function NedelecRefFE(::Type{et},p::Polytope,order::Integer) where et
  D = num_dims(p)

  if is_n_cube(p)
    prebasis = QGradMonomialBasis{D}(et,order) # Prebasis
    eb = MonomialBasis{1}(et,order)            # Edge basis
    fb = QGradMonomialBasis{D-1}(et,order-1)   # Face basis
    cb = QCurlGradMonomialBasis{D}(et,order-1) # Cell basis
  elseif is_simplex(p)
    prebasis = Polynomials.NedelecPrebasisOnSimplex{D}(order) # Prebasis
    eb = MonomialBasis{1}(et,order)                           # Edge basis
    fb = MonomialBasis{D-1}(VectorValue{D-1,et},order-1,Polynomials._p_filter) # Face basis
    cb = MonomialBasis{D}(VectorValue{D,et},order-D+1,Polynomials._p_filter)   # Cell basis
  else
    @unreachable "Nedelec Reference FE only implemented for n-cubes and simplices"
  end

  function cmom(φ,μ,ds) # Cell moment function: σ_K(φ,μ) = ∫(φ⋅μ)dK
    Broadcasting(Operation(⋅))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function: σ_F(φ,μ) = ∫((φ×n)⋅μ)dF
    o = get_facet_orientations(ds.poly)[ds.face] # This is a hack to avoid a sign map
    n = o*get_facet_normal(ds)
    E = get_extension(ds)
    Eμ = Broadcasting(Operation(⋅))(E,μ) # We have to extend the basis to 3D
    φn = Broadcasting(Operation(×))(n,φ)
    Broadcasting(Operation(⋅))(φn,Eμ)
  end
  function emom(φ,μ,ds) # Edge moment function: σ_E(φ,μ) = ∫((φ⋅t)*μ)dE
    t = get_edge_tangent(ds)
    φt = Broadcasting(Operation(⋅))(φ,t)
    Broadcasting(Operation(*))(φt,μ)
  end

  if D == 2
    moments = [ # In 2D we do not have face moments
      (get_dimrange(p,1),emom,eb),   # Edge moments
      (get_dimrange(p,D),cmom,cb)    # Cell moments
    ]
  else
    @assert D == 3
    moments = [
      (get_dimrange(p,1),emom,eb),   # Edge moments
      (get_dimrange(p,D-1),fmom,fb), # Face moments
      (get_dimrange(p,D),cmom,cb)    # Cell moments
    ]
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

struct CoVariantPiolaMap <: Map end

function evaluate!(
  cache,
  ::Broadcasting{typeof(∇)},
  a::Fields.BroadcastOpFieldArray{CoVariantPiolaMap})
  v, Jt = a.args
  # Assuming J comes from an affine map
  ∇v = Broadcasting(∇)(v)
  k = CoVariantPiolaMap()
  Broadcasting(Operation(k))(∇v,Jt)
end

function lazy_map(
  ::Broadcasting{typeof(gradient)},
  a::LazyArray{<:Fill{Broadcasting{Operation{CoVariantPiolaMap}}}})
  v, Jt = a.args
  ∇v = lazy_map(Broadcasting(∇),v)
  k = CoVariantPiolaMap()
  lazy_map(Broadcasting(Operation(k)),∇v,Jt)
end

function evaluate!(cache,::CoVariantPiolaMap,v::Number,Jt::Number)
  v⋅transpose(inv(Jt)) # we multiply by the right side to compute the gradient correctly
end

function evaluate!(cache,k::CoVariantPiolaMap,v::AbstractVector{<:Field},phi::Field)
  Jt = ∇(phi)
  Broadcasting(Operation(k))(v,Jt)
end

function lazy_map(
  k::CoVariantPiolaMap,
  cell_ref_shapefuns::AbstractArray{<:AbstractArray{<:Field}},
  cell_map::AbstractArray{<:Field})

  cell_Jt = lazy_map(∇,cell_map)
  lazy_map(Broadcasting(Operation(k)),cell_ref_shapefuns,cell_Jt)
end
