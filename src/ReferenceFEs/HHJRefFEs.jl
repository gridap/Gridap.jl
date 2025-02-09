
struct HellanHerrmannJhonson <: PushforwardRefFE end

const hhj = HellanHerrmannJhonson()

Pushforward(::Type{<:HellanHerrmannJhonson}) = DoubleContraVariantPiolaMap()

"""
    struct HellanHerrmannJhonson <: PushforwardRefFE end
    HellanHerrmannJhonsonRefFE(::Type{T},p::Polytope,order::Integer) where T

Hellan-Herrmann-Jhonson reference finite element.

References:

- `The Hellan-Herrmann-Johnson method with curved elements`, Arnold and Walker (2020)

"""
function HellanHerrmannJhonsonRefFE(::Type{T},p::Polytope,order::Integer) where T
  @assert p == TRI "HellanHerrmannJhonson Reference FE only defined for TRIangles"
  
  VT = SymTensorValue{2,T}
  prebasis = MonomialBasis(Val(2),VT,order,Polynomials._p_filter)
  #fb = MonomialBasis(Val(1),T,order,Polynomials._p_filter)
  fb = get_shapefuns(LagrangianRefFE(T, SEGMENT, order))
  #cb = MonomialBasis(Val(2),VT,order-1,Polynomials._p_filter)
  #cb = get_shapefuns(LagrangianRefFE(VT, TRI, order-1))
  cb = map(constant_field,[VT(0.,1.,0.),VT(-2.,1.,0.),VT(0.,-1.,2.)])

  function cmom(φ,μ,ds) # Cell and Node moment function: σ_K(φ,μ) = ∫(φ:μ)dK
    Broadcasting(Operation(⊙))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function (normal) : σ_F(φ,μ) = ∫((n·φ·n)*μ)dF
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    nφn = Broadcasting(Operation(⋅))(n,φn)
    Broadcasting(Operation(*))(nφn,μ)
  end

  moments = Tuple[
    (get_dimrange(p,1),fmom,fb), # Face moments
  ]
  if order > 0
    push!(moments,(get_dimrange(p,2),cmom,cb))  # Cell moments
  end

  return MomentBasedReferenceFE(HellanHerrmannJhonson(),p,prebasis,moments,DivConformity())
end

Polynomials.get_order(f::Fields.LinearCombinationFieldVector) = get_order(f.fields)
Polynomials.get_order(f::AbstractVector{<:ConstantField}) = 0

function ReferenceFE(p::Polytope,::HellanHerrmannJhonson, order)
  HellanHerrmannJhonsonRefFE(Float64,p,order)
end

function ReferenceFE(p::Polytope,::HellanHerrmannJhonson,::Type{T}, order) where T
  HellanHerrmannJhonsonRefFE(T,p,order)
end

function Conformity(reffe::GenericRefFE{HellanHerrmannJhonson},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a HellanHerrmannJhonson reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
      """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{HellanHerrmannJhonson}, conf::DivConformity)
  reffe.face_dofs
end

function get_face_dofs(reffe::GenericRefFE{HellanHerrmannJhonson,Dc}) where Dc
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
