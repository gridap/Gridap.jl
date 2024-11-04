
using Gridap
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.Arrays
using Gridap.Geometry

using FillArrays

using Gridap.ReferenceFEs: MonomialBasis, JacobiBasis
using Gridap.Polynomials, Gridap.TensorValues

############################################################################################
# Fixes

function Arrays.return_cache(k::Fields.BroadcastOpFieldArray{typeof(∘)},x::AbstractArray{<:Point})
  f, g = k.args
  cg = return_cache(g,x)
  gx = evaluate!(cg,g,x)
  cf = return_cache(f,gx)
  return cg, cf
end

function Arrays.evaluate!(cache, k::Fields.BroadcastOpFieldArray{typeof(∘)}, x::AbstractArray{<:Point})
  cg, cf = cache
  f, g = k.args
  gx = evaluate!(cg,g,x)
  fgx = evaluate!(cf,f,gx)
  return fgx
end

function Arrays.return_value(k::Fields.OperationField{typeof(∘)},x::AbstractArray{<:Point})
  f, g = k.fields
  gx  = return_value(g,x)
  fgx = return_value(f,gx)
  return fgx
end

function Arrays.return_cache(k::Fields.OperationField{typeof(∘)},x::AbstractArray{<:Point})
  f, g = k.fields
  cg = return_cache(g,x)
  gx = evaluate!(cg,g,x)
  cf = return_cache(f,gx)
  return cg, cf
end

function Arrays.evaluate!(cache, k::Fields.OperationField{typeof(∘)}, x::AbstractArray{<:Point})
  cg, cf = cache
  f, g = k.fields
  gx = evaluate!(cg,g,x)
  fgx = evaluate!(cf,f,gx)
  return fgx
end

Arrays.return_type(::Polynomials.QGradMonomialBasis{D,T}) where {D,T} = T

############################################################################################
# FaceMeasure

struct FaceMeasure{Df,Dc}
  poly::Polytope{Dc}
  face::Int
  quad::Quadrature
  function FaceMeasure{Df}(poly::Polytope{Dc},face::Int,order::Int) where {Df,Dc}
    fpoly = get_face_polytope(poly,Df,face)
    quad = Quadrature(fpoly,order)
    new{Df,Dc}(poly,face,quad)
  end
end

function get_face_polytope(p::Polytope{Dc},Df::Int,face::Int) where Dc
  return get_reffaces(p)[get_face_type(p)[get_offset(p,Df) + face]]
end

function get_normal(m::FaceMeasure{Df,Dc}) where {Df,Dc}
  @assert Df == Dc - 1
  n = get_facet_normal(m.poly)
  return ConstantField(n[m.face])
end

function get_tangent(m::FaceMeasure{1,Dc}) where {Dc}
  t = get_edge_tangent(m.poly)
  return ConstantField(t[m.face])
end

function Geometry.get_cell_map(m::FaceMeasure{Df,Dc}) where {Df,Dc}
  if Df == Dc
    return GenericField(identity)
  end
  fp = get_face_polytope(m.poly,Df,m.face)
  coords = get_face_coordinates(m.poly,Df)[m.face]
  fields = get_shapefuns(LagrangianRefFE(Float64,fp,1))
  return linear_combination(coords,fields)
end

function get_extension(m::FaceMeasure{Df,Dc}) where {Df,Dc}
  @assert Df == Dc - 1
  vs = ReferenceFEs._nfaces_vertices(Float64,m.poly,Df)[m.face]
  return ConstantField(TensorValue(hcat([vs[2]-vs[1]...],[vs[3]-vs[1]...])))
end

function Arrays.evaluate(f,ds::FaceMeasure)
  x = get_coordinates(ds.quad)
  w = get_weights(ds.quad)
  fx = evaluate(f,x)
  return w .* fx, x
end

############################################################################################
# Constant basis: Basis for a tensor type

constant_basis(T::Type{<:Real}) = [one(T)]
function constant_basis(V::Type{<:MultiValue})
  T = eltype(V)
  n = num_components(V)
  z, o = zero(T), one(T)
  return [V(ntuple(i -> ifelse(i == j, o, z),Val(n))) for j in 1:n]
end

############################################################################################
using Gridap.ReferenceFEs: ReferenceFEName, Conformity
struct MomentBasedReffe{T<:ReferenceFEName} <: ReferenceFEName
  name :: T
end

"""
A moment is given by a triplet (f,σ,μ) where 
  - f is id of a face Fk
  - σ is a function σ(φ,μ,ds) that returns a Field-like object to be integrated over Fk
  - μ is a polynomials basis on Fk

Open questions: 
  - Do we want to keep having structures face -> data? I guess if we had more than a single 
    moment per face, we would aggregate them. 
  - Can we always determine the minimum integration order for each moment? 
"""
function MomentBasedReferenceFE(
  name::ReferenceFEName,
  p::Polytope{D},
  prebasis::AbstractVector{<:Field},
  moments::AbstractVector{<:Tuple{<:Integer,<:Function,<:AbstractArray{<:Field}}},
  conformity::Conformity;
) where D

  # TODO:  Basis of the constants for the tensor-type we have
  T = return_type(prebasis)
  order = get_order(prebasis)
  φ = constant_basis(VectorValue{D,T})

  # TODO: This has to be something that can fully contract with the prebasis φ
  V = VectorValue{D,Float64}

  # TODO: Do we want these of length n_moments or n_faces?
  n_faces = num_faces(p)
  n_moments = length(moments)
  face_moments = Vector{Array{V}}(undef, n_moments)
  face_nodes = Vector{UnitRange{Int}}(undef, n_moments)
  face_own_dofs = [Int[] for _ in 1:n_faces]
  nodes = Point{D}[]

  k = 1
  n_nodes = 1
  n_dofs = 1
  for (face,σ,μ) in moments
    d = get_facedims(p)[face]
    lface = face - get_offset(p,d)

    qdegree = order + get_order(μ) + 1
    ds = FaceMeasure{d}(p, lface, qdegree)

    fmap = get_cell_map(ds)
    φf = transpose(Broadcasting(Operation(∘))(map(constant_field,φ),fmap))
    vals, f_coords = evaluate(σ(φf,μ,ds),ds)
    coords = evaluate(fmap,f_coords)

    face_moments[k] = map(v -> v⋅φ, eachslice(vals, dims=(1,2))) # (nN, nμ, nφ) ⋅ nφ -> (nN, nμ)
    face_nodes[k] = n_nodes:(n_nodes+length(coords)-1)
    append!(nodes, coords)
    append!(face_own_dofs[face], n_dofs:(n_dofs+size(vals,2)-1))

    k += 1
    n_nodes += length(coords)
    n_dofs += size(vals,2)
  end

  dof_basis = MomentBasedDofBasis(nodes, face_moments, face_nodes)
  metadata = nothing

  return GenericRefFE{typeof(MomentBasedReffe(name))}(
    n_dofs,
    p,
    prebasis,
    dof_basis,
    conformity,
    metadata,
    face_own_dofs
  )
end

function cmom(φ,μ,ds)
  Broadcasting(Operation(⋅))(φ,μ)
end

function fmom_dot(φ,μ,ds)
  n = get_normal(ds)
  φn = Broadcasting(Operation(⋅))(φ,n)
  Broadcasting(Operation(*))(φn,μ)
end

function fmom_cross(φ,μ,ds)
  n = get_normal(ds)
  E = get_extension(ds)
  Eμ = Broadcasting(Operation(⋅))(E,μ) # We have to extend the basis to 3D (see Nedelec)
  φn = Broadcasting(Operation(×))(φ,n)
  Broadcasting(Operation(⋅))(φn,Eμ)
end

function emom(φ,μ,ds)
  t = get_tangent(ds)
  φt = Broadcasting(Operation(⋅))(φ,t)
  Broadcasting(Operation(*))(φt,μ)
end

# RT implementation

D = 2
p = (D==2) ? QUAD : HEX
order = 1

prebasis = QCurlGradMonomialBasis{D}(Float64,order)
cb = QGradJacobiPolynomialBasis{D}(Float64,order-1)
fb = JacobiBasis(Float64,SEGMENT,order)
moments = [
  [(f+get_offset(p,1),fmom_dot,fb) for f in 1:num_faces(p,1)]..., # Face moments
  (num_faces(p),cmom,cb) # Cell moments
]
reffe = MomentBasedReferenceFE(RaviartThomas(),p,prebasis,moments,DivConformity())
dofs = get_dof_basis(reffe)

rt_reffe = RaviartThomasRefFE(Float64,p,order)
rt_dofs = get_dof_basis(rt_reffe)

Mrt = evaluate(rt_dofs,prebasis)
M = evaluate(dofs,prebasis)
M == Mrt

# ND implementation

D = 2
p = (D==2) ? QUAD : HEX
order = 1

prebasis = QGradMonomialBasis{D}(Float64,order)
cb = QCurlGradMonomialBasis{D}(Float64,order-1)
fb = QGradMonomialBasis{D-1}(Float64,order-1)
eb = MonomialBasis(Float64,SEGMENT,order)
moments = [
  [(f+get_offset(p,1),emom,eb) for f in 1:num_faces(p,1)]..., # Edge moments
  (num_faces(p),cmom,cb) # Cell moments
]
reffe = MomentBasedReferenceFE(Nedelec(),p,prebasis,moments,CurlConformity())
dofs = get_dof_basis(reffe)

nd_reffe = NedelecRefFE(Float64,p,order)
nd_dofs = get_dof_basis(nd_reffe)

Mnd = evaluate(nd_dofs,prebasis)
M = evaluate(dofs,prebasis)
M == Mnd

# 3D ND implementation

D = 3
p = (D==2) ? QUAD : HEX
order = 1

prebasis = QGradMonomialBasis{D}(Float64,order)
cb = QCurlGradMonomialBasis{D}(Float64,order-1)
fb = QGradMonomialBasis{D-1}(Float64,order-1)
eb = MonomialBasis(Float64,SEGMENT,order)
moments = [
  [(f+get_offset(p,1),emom,eb) for f in 1:num_faces(p,1)]..., # Edge moments
  [(f+get_offset(p,2),fmom_cross,fb) for f in 1:num_faces(p,2)]..., # Face moments
  (num_faces(p),cmom,cb) # Cell moments
]
reffe = MomentBasedReferenceFE(Nedelec(),p,prebasis,moments,CurlConformity())
dofs = get_dof_basis(reffe)

nd_reffe = NedelecRefFE(Float64,p,order)
nd_dofs = get_dof_basis(nd_reffe)

Mnd = evaluate(nd_dofs,prebasis)
M = evaluate(dofs,prebasis)
M == Mnd # This is because the get_facet_orientations hack in Nedelec... Why is it necessary? 

############################################################################################

#φ = MonomialBasis(VectorValue{2,Float64},QUAD,1)
φ = ConstantField(VectorValue(ntuple(i -> 1.0, 2)))

dc = FaceMeasure{2}(QUAD,1,2)
mc = get_cell_map(dc)
φc = transpose(Broadcasting(Operation(∘))(φ,mc))
μc = MonomialBasis(VectorValue{2,Float64},QUAD,1)

df = FaceMeasure{1}(QUAD,1,2)
mf = get_cell_map(df)
φf = transpose(Broadcasting(Operation(∘))(φ,mf))
μf = MonomialBasis(Float64,SEGMENT,1)

σc = cmom(φc,μc,dc)
σf = fmom(φf,μf,df)

vc = evaluate(σc,dc)
vf = evaluate(σf,df)

sum(vc, dims=1)



"""
List of issues: 

- ReferenceFEs is loaded before CellData, i.e we do NOT have access to the 
  CellField machinery to compute the moments.

- Most operations that are defined for CellFields are not 100% working for arrays of Fields, 
  where we tend to use the Broadcasting + Operation machinery.
  For example, ∇(φ) is explicitly deactivated in favor of Broadcasting(∇)(φ).
"""

p = QUAD

pts = [Point(0.0,0.0), Point(1.0,0.0), Point(1.0,1.0), Point(0.0,1.0)]

φ = ReferenceFEs.MonomialBasis(Float64,p,1)
Broadcasting(∇)(φ)


evaluate(*,φ,φ)
evaluate(+,φ,φ)
evaluate(-,φ,φ)

evaluate(evaluate(Broadcasting(+),φ,φ),pts)
evaluate(Broadcasting(Operation(+))(φ,φ),pts)

evaluate(evaluate(⋅,φ,φ),pts)
evaluate(evaluate(Broadcasting(Operation(⋅)),φ,φ),pts)


############################################################################################
# Deprecated code

function MomentBasedReferenceFE(
  p::Polytope{D},
  order::Int,
  face_moments::AbstractVector{<:Union{Function,Nothing}},
  face_basis::AbstractVector{<:Union{AbstractVector{<:Field},Nothing}}
) where D
  _face_moments = Vector{AbstractVector{<:Union{Function,Nothing}}}(undef, D+1)
  _face_basis = Vector{AbstractVector{<:Union{AbstractVector{<:Field},Nothing}}}(undef, D+1)
  for d in 0:D
    _face_moments[d+1] = Fill(face_moments[d+1], num_faces(p,d))
    _face_basis[d+1] = Fill(face_basis[d+1], num_faces(p,d))
  end
  MomentBasedReferenceFE(
    p,
    order,
    _face_moments,
    _face_basis
    #[Fill(m,num_faces(p,d-1)) for (d,m) in enumerate(face_moments)],
    #[Fill(b,num_faces(p,d-1)) for (d,b) in enumerate(face_basis)],
  )
end

function MomentBasedReferenceFE(
  p::Polytope{D},
  order::Int,
  face_moments::AbstractVector{<:AbstractVector{<:Union{Function,Nothing}}},
  face_basis::AbstractVector{<:AbstractVector{<:Union{AbstractVector{<:Field},Nothing}}}
) where D
  @assert length(face_moments) == length(face_basis) == D+1
  @assert all(length(face_moments[d+1]) == length(face_basis[d+1]) for d in 0:D)

  d_to_face_to_coords = [Vector{Vector{Point{D}}}(undef, num_faces(p,d)) for d in 0:D]
  d_to_face_to_vals = [Vector{Matrix{Float64}}(undef, num_faces(p,d)) for d in 0:D]
  for d in 0:D
    for face in 1:num_faces(p,d)
      σ = face_moments[d+1][face]
      μ = face_basis[d+1][face]

      if !isnothing(σ)
        @assert !isnothing(μ)
        qdegree = order*get_order(μ)
        ds = FaceMeasure{d}(p, face, qdegree)

        fmap = get_cell_map(ds)
        φ = Operation(∘)(ConstantField(VectorValue(ntuple(i -> 1.0, D))),fmap)
        vals, coords = evaluate(σ(φ,μ,ds),ds)

        d_to_face_to_coords[d+1][face] = evaluate(fmap,coords)
        d_to_face_to_vals[d+1][face] = vals
      end
    end
  end

  return d_to_face_to_coords, d_to_face_to_vals
end

