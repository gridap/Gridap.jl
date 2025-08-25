
module GeometricDecompositionTest

using Test
using Gridap.Helpers
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Fields: MockField
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.Arrays
using Gridap.Visualization

using StaticArrays
using LinearAlgebra

using Gridap
using Gridap.FESpaces

using Gridap.ReferenceFEs: FaceMeasure, set_face!

import Gridap.Arrays: return_cache
import Gridap.Arrays: evaluate!

##########################
# FaceIntegralFormVector #
##########################

# The MomentBasedDofBasis cannot implement the L2 trace norms on
# Polytope faces, because the norm is nonlinear. An alternative implementation is
# introduced here to test the trace norms of the geometrically decomposed bases.

const _μ = MonomialBasis(Val(0),Float64,0) # a length-1 dummy basis that won't be evaluated

"""
    abstract type Form <: Map

Abstract type for a form over a functional space (typically a polynomial space).
The domain is a [`Field`](@ref) set and the range the scalar set.

If put in src, this should be a supertype of `Dof`, the latter being a *linear* `Form`.
"""
abstract type Form <: Map end

"""
    FaceIntegralFormVector(
      p::Polytope{D},
      order::Int,
      faces_and_forms_integrands::Vector{Tuple{Vector{Int},Function}}
    )

Constructs a vector of (nonlinear) forms that are integrals over faces of `p`.

`faces_and_forms_integrands` is a collection of form definitions, each one is
given by a couple (faces_ids, σ) where
  - faces_ids is vector of ids of faces Fₖ of `p`
  - σ is the form integrand φ,_,ds -> σ(φ,_,ds) that returns a scalar-valued Field-like object to be integrated over each Fₖ

The forms are thus defined by φ -> ∫_Fₖ σ(φ,_,ds).

In the final basis, the forms are ordered by moment descriptor, then by face.

We are assuming that all the faces in a moment are of the same type.
"""
struct FaceIntegralFormVector{T,FI,DS} <: AbstractVector{Form}
  form_faces::Vector{Vector{NTuple{2,Int}}} # f -> faces over which f is integrated computed
  form_integrands::FI # i -> σᵢ
  form_measures::DS # tuple of `FaceMeasure`s
  face_own_forms::Vector{Vector{Int}} # "inverse" of form_faces

  function FaceIntegralFormVector(faces,integrands,measures,f_own_forms)
    T = Float64 # TODO
    FI = typeof(integrands)
    DS = typeof(measures)
    new{T,FI,DS}(faces,integrands,measures,f_own_forms)
  end

  function FaceIntegralFormVector(p::Polytope,order::Int,forms)
    n_faces = num_faces(p)
    n_form = length(forms)
    face_dims = get_facedims(p)
    face_offsets = get_offsets(p)
    reffaces, face_types = ReferenceFEs._compute_reffaces_and_face_types(p)

    # Create facemeasures and integrand for each integral form
    # Count number of forms per face
    # compute local to global form index
    integrands = ()
    measures = ()
    face_n_forms = zeros(Int,n_faces)
    face_form_ids = Vector{Vector{NTuple{2,Int}}}(undef,n_form)
    form_I = 0
    for (k,(faces,σ)) in enumerate(forms)
      face_ids = Vector{NTuple{2,Int}}(undef, length(faces))
      for (i,face) in enumerate(faces)
        d = face_dims[face]
        lface = face - face_offsets[d+1]
        face_ids[i] = (form_I + i, lface)
      end
      face_form_ids[k] = face_ids
      form_I += length(faces)

      face_n_forms[faces] .+= 1
      integrands = (integrands..., σ)

      ftype = face_types[first(faces)]
      @check all(isequal(ftype), face_types[faces])
      fp = reffaces[ftype]
      measures = (measures..., FaceMeasure(p,fp,order) )
    end

    # Compute face forms and nodes indices
    n_forms = 0
    face_own_forms = Vector{Vector{Int}}(undef,n_faces)
    for face in 1:n_faces
      n_forms_i = face_n_forms[face]
      face_own_forms[face] = collect((n_forms+1):(n_forms+n_forms_i))
      n_forms += n_forms_i
    end

    FaceIntegralFormVector(face_form_ids,integrands,measures,face_own_forms)
  end
end

Base.size(a::FaceIntegralFormVector) = (num_forms(a),)
Base.axes(a::FaceIntegralFormVector) = (Base.OneTo(num_forms(a)),)
Base.getindex(::FaceIntegralFormVector,::Integer) = Form()
Base.IndexStyle(::FaceIntegralFormVector) = IndexLinear()
num_forms(b::FaceIntegralFormVector) = mapreduce(length, +, b.face_own_forms)

function return_cache(b::FaceIntegralFormVector{T}, f) where {T}
  cms = ()
  for (_,σ,ds) in zip(b.form_faces, b.form_integrands, b.form_measures)
    cms = ( cms..., return_cache(σ,f,_μ,ds))
  end

  r = Array{T}(undef, (num_forms(b), size(f)...))
  c = CachedArray(r)
  return cms, c
end

function evaluate!(cache, b::FaceIntegralFormVector, f::AbstractVector{<:Field})
  cms, c = cache
  nf = length(f)
  setsize!(c, (length(b), nf))

  forms_vals = c.array
  for (cds,faces,σ,ds) in zip(cms, b.form_faces, b.form_integrands, b.form_measures)
    for (form_id, lface) in faces
      set_face!(ds,lface)
      # vals: v_{i,1,j} =  σ(f_j(n_i), _μ₁, ds), size (nN, 1, nf)
      vals, _ = evaluate!(cds,σ,f,_μ,ds)
      forms_vals[form_id,:] = reshape(sum(vals,dims=1), (nf,))
    end
  end
  forms_vals
end


###########################################################################################
# Traces L2-norm descriptors and geometric decompositions trace-forms for each conformity #
###########################################################################################

function zero_moment(φ,_,ds) # zero_f(φ) = ∫ 0 df
  C0 = ConstantField(0)
  φ0 = Broadcasting(Operation(*))(C0,φ)
  Broadcasting(Operation(norm))(φ0)
end

function HGrad_face_tr_norm(φ,_,ds) # tr_f(φ) = ∫ ‖φ‖ df
  Broadcasting(Operation(norm))(φ)
end

function HCurl_edge_tr_norm(φ,_,ds) # tr_E(φ) = ∫ ‖φ·t‖ dE
  t = get_edge_tangent(ds)
  φt = Broadcasting(Operation(⋅))(t,φ)
  Broadcasting(Operation(norm))(φt)
end

function HCurl_facet_tr_norm(φ,_,ds) # tr_F(φ) = ∫ ‖φ×n‖ dF
  n = get_facet_normal(ds)
  φt = Broadcasting(Operation(×))(φ,n)
  Broadcasting(Operation(norm))(φt)
end

function HDiv_facet_tr_norm(φ,_,ds) # tr_F(φ) = ∫ ‖φ⋅n‖ dF
  n = get_facet_normal(ds)
  φn = Broadcasting(Operation(⋅))(φ,n)
  Broadcasting(Operation(norm))(φn)
end


# zero form on each face
function get_trace_forms(p::Polytope{D}, field, ::L2Conformity) where D
  face_ranges = get_dimranges(p)
  trace_forms = Tuple[ (f_range, zero_moment) for f_range in face_ranges ]
  qorder = 2*get_order(field)+1
  FaceIntegralFormVector(p,qorder,trace_forms)
end

# norm form on each boundary face, zero inside
function get_trace_forms(p::Polytope{D}, field, ::GradConformity) where D
  face_ranges = get_dimranges(p)
  trace_forms = Tuple[ (face_ranges[d], HGrad_face_tr_norm) for d in 1:D ]
  push!(trace_forms, (last(face_ranges), zero_moment))
  qorder = 2*get_order(field)+1
  FaceIntegralFormVector(p,qorder,trace_forms)
end

# zero form on vertices and inside, norm of tangential trace on edges (φ⋅t) and facets (φ×n)
function get_trace_forms(p::Polytope{D}, field, ::CurlConformity) where D
  D>3 && @notimplemented "tangential trace to 2D faces in D>3 notimplemented"
  face_ranges = get_dimranges(p)
  trace_forms = Tuple[ (first(face_ranges), zero_moment) ]
  D≥2 && push!(trace_forms, (face_ranges[2], HCurl_edge_tr_norm))
  D≥3 && push!(trace_forms, (face_ranges[D], HCurl_facet_tr_norm))
  push!(  trace_forms, (last( face_ranges), zero_moment))
  qorder = 2*get_order(field)+1
  FaceIntegralFormVector(p,qorder,trace_forms)
end

# zero form on vertices, edges and inside, norm of normal trace to facets (φ⋅n)
function get_trace_forms(p::Polytope{D}, field, ::DivConformity) where D
  face_ranges = get_dimranges(p)
  trace_forms = Tuple[ (face_ranges[d], zero_moment) for d in 1:D-1 ]
  push!(trace_forms, (face_ranges[D], HDiv_facet_tr_norm))
  push!(trace_forms, (last( face_ranges), zero_moment))
  qorder = 2*get_order(field)+1
  FaceIntegralFormVector(p,qorder,trace_forms)
end


# Traces Test function
function _test_geometric_decomposition(b,p,conf,
  face_own_funs=get_face_own_funs(b,p,conf), skip_check=false)

  if !skip_check
    @test has_geometric_decomposition(b,p,conf)
  end

  faces = get_faces(p)
  tr_forms = get_trace_forms(p,b,conf) # one form for each face in faces
  b_face_traces = evaluate(tr_forms, b)
  tr_iszero = map(b_tr -> abs(b_tr)<tol, b_face_traces)

  pass = true
  for (f, own_funs) in enumerate(face_own_funs)
    f_is_in_g = broadcast(∈, f, faces)
    for fun in own_funs
      tr_fun_iszero_on_g = tr_iszero[:,fun]
      pass = pass && all(tr_fun_iszero_on_g .|| f_is_in_g)

      #if !pass # for debug
      #  @show b, p
      #  @show faces, face_own_funs
      #  @show "function ", fun, " face ", f
      #  @show b_face_traces[:,fun], "trace"
      #  @show tr_fun_iszero_on_g, "is_zero"
      #  @show f_is_in_g, "is f in g "
      #  @unreachable
      #end
    end
  end
  @test pass
end

r = 3
et = Float64
tol = 1000*eps(et)

SIMPL4 = ExtrusionPolytope(tfill(TET_AXIS,Val(4)))
NCUBE4 = ExtrusionPolytope(tfill(HEX_AXIS,Val(4)))

##################################
# Grad conforming decompositions #
##################################

conf = GradConformity()

for p in (SEGMENT,TRI,TET,SIMPL4)
  D = num_dims(p)

  b = BernsteinBasisOnSimplex(Val(D),et,r)
  _test_geometric_decomposition(b,p,conf)

  b = BernsteinBasisOnSimplex(Val(D),SkewSymTensorValue{3,et},r)
  _test_geometric_decomposition(b,p,conf)

  b = BarycentricPmΛBasis(Val(D),et,r,0)
  _test_geometric_decomposition(b,p,conf)

  b = BarycentricPΛBasis(Val(D),et,r,0)
  _test_geometric_decomposition(b,p,conf)
end

ModalC0 = Polynomials.ModalC0

# non-zero trace on both x=0. and 1. sides
b = CartProdPolyBasis(ModalC0,Val(2),et,(0,1))
@test !has_geometric_decomposition(b,QUAD,conf)

for p in (SEGMENT,QUAD,HEX,NCUBE4)
  D = num_dims(p)

  b = CartProdPolyBasis(ModalC0,  Val(D),et,r)
  _test_geometric_decomposition(b,p,conf)

  b = CartProdPolyBasis(Bernstein,Val(D),et,r)
  _test_geometric_decomposition(b,p,conf)

  b = CartProdPolyBasis(ModalC0,  Val(D), SkewSymTensorValue{3,et}, r)
  _test_geometric_decomposition(b,p,conf)

  b = CartProdPolyBasis(Bernstein,Val(D), SkewSymTensorValue{3,et}, r)
  _test_geometric_decomposition(b,p,conf)
end


#########################################
# Curl Conform Geometric decompositions #
#########################################

conf = CurlConformity()

for p in (SEGMENT,TRI,TET) #,SIMPL4)
  D = num_dims(p)
  k = 1

  b = BarycentricPmΛBasis(Val(D),et,r,k) # Nédélec
  _test_geometric_decomposition(b,p,conf)

  b = BarycentricPΛBasis(Val(D),et,r,k)  # Nédélec 2nd kind
  _test_geometric_decomposition(b,p,conf)
end

D=2
V = VectorValue{D,et}
orders = [ (i==j ? 0 : 1) for i in 1:D, j in 1:D ]
b = CompWiseTensorPolyBasis{D}(ModalC0,V,orders)   # RT
@test !has_geometric_decomposition(b,QUAD,GradConformity())
@test !has_geometric_decomposition(b,QUAD,DivConformity())

for p in (QUAD,HEX)
  D = num_dims(p)
  V = VectorValue{D,et}
  orders = [ r-1 + (i==j ? 0 : 1) for i in 1:D, j in 1:D ]

  b = CompWiseTensorPolyBasis{D}(ModalC0,V,orders)   # Nédélec
  _test_geometric_decomposition(b,p,conf)

  b = CompWiseTensorPolyBasis{D}(Bernstein,V,orders) # Nédélec
  _test_geometric_decomposition(b,p,conf)
end

########################################
# Div Conform Geometric decompositions #
########################################

conf = DivConformity()

for p in (TRI,TET,SIMPL4)
  D = num_dims(p)
  k = D-1 # flux forms
  rotate_90 = D==2

  b = BarycentricPmΛBasis(Val(D),et,r,k; rotate_90) # RT
  _test_geometric_decomposition(b,p,conf)

  b = BarycentricPΛBasis(Val(D),et,r,k; rotate_90)  # BDM
  _test_geometric_decomposition(b,p,conf)
end

D=2
V = VectorValue{D,et}
orders = [ (i==j ? 1 : 0) for i in 1:D, j in 1:D ]
b = CompWiseTensorPolyBasis{D}(ModalC0,V,orders)   # RT
@test !has_geometric_decomposition(b,QUAD,GradConformity())
@test !has_geometric_decomposition(b,QUAD,CurlConformity())

for p in (QUAD,HEX)
  D = num_dims(p)
  V = VectorValue{D,et}
  orders = [ r-1 + (i==j ? 1 : 0) for i in 1:D, j in 1:D ]

  b = CompWiseTensorPolyBasis{D}(ModalC0,V,orders)   # RT
  _test_geometric_decomposition(b,p,conf)

  b = CompWiseTensorPolyBasis{D}(Bernstein,V,orders) # RT
  _test_geometric_decomposition(b,p,conf)
end

#######################################
# L2 Conform Geometric decompositions #
#######################################

conf = L2Conformity()

p = VERTEX
b = [ ConstantField(1) ]
@test has_geometric_decomposition(b,p,conf)
@test get_face_own_funs(b,p,conf) == [ [1] ]

p = QUAD
b = [ ConstantField(vi) for vi in component_basis(VectorValue{3,et}) ]
@test has_geometric_decomposition(b,p,conf)
@test get_face_own_funs(b,p,conf) == [
  Int[], Int[], Int[], Int[], Int[], Int[], Int[], Int[], Int[1,2,3],
]

# Other APIs
conf = CDConformity( (CONT,DISC) )
@test !has_geometric_decomposition(b,p,conf)
@test_throws ErrorException get_face_own_funs(b,p,conf)

end # module

