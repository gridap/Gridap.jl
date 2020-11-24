
abstract type Conformity end

struct L2Conformity <: Conformity end

"""
    abstract type ReferenceFE{D} <: GridapType

Abstract type representing a Reference finite element. `D` is the underlying coordinate space dimension.
We follow the Ciarlet definition. A reference finite element
is defined by a polytope (cell topology), a basis of an interpolation space
of top of this polytope (denoted here as the prebasis), and a basis of the dual of this space
(i.e. the degrees of freedom). From this information one can compute the shape functions
(i.e, the canonical basis of w.r.t. the degrees of freedom) with a simple change of basis.
In addition, we also encode in this type information about how the interpolation space
in a reference finite element is "glued" with neighbors in order to build conforming
cell-wise spaces.

The `ReferenceFE` interface is defined by overloading these methods:

- [`num_dofs(reffe::ReferenceFE)`](@ref)
- [`get_polytope(reffe::ReferenceFE)`](@ref)
- [`get_prebasis(reffe::ReferenceFE)`](@ref)
- [`get_dof_basis(reffe::ReferenceFE)`](@ref)
- [`Conformity(reffe::ReferenceFE)`](@ref)
- [`get_face_own_dofs(reffe::ReferenceFE,conf::Conformity)`](@ref)
- [`get_face_own_dofs_permutations(reffe::ReferenceFE,conf::Conformity)`](@ref)
- [`get_face_dofs(reffe::ReferenceFE)`](@ref)

The interface is tested with
- [`test_reference_fe(reffe::ReferenceFE)`](@ref)

"""
abstract type ReferenceFE{D} <: GridapType end

# Extensible factory function
function ReferenceFE(p::Polytope,basis::Symbol,args...;kwargs...)
  ReferenceFE(p,Val(basis),args...;kwargs...)
end

function  ReferenceFE(p::Polytope,::Val{T},args...;kwargs...) where T
  @unreachable """\n
  Undefined factory function ReferenceFE for symbol $T
  """
end

ReferenceFE(basis::Symbol,args...;kwargs...) = (basis, args, kwargs)

"""
    num_dofs(reffe::ReferenceFE) -> Int

Returns the number of DOFs.
"""
function num_dofs(reffe::ReferenceFE)
  @abstractmethod
end

"""
    get_polytope(reffe::ReferenceFE) -> Polytope

Returns the underlying polytope object.
"""
function get_polytope(reffe::ReferenceFE)
  @abstractmethod
end

"""
    get_prebasis(reffe::ReferenceFE) -> Field

Returns the underlying prebasis encoded as a `Field` object.
"""
function get_prebasis(reffe::ReferenceFE)
  @abstractmethod
end

get_order(reffe::ReferenceFE) = get_order(get_prebasis(reffe))

"""
    get_dof_basis(reffe::ReferenceFE) -> Dof

Returns the underlying dof basis encoded in a `Dof` object.
"""
function get_dof_basis(reffe::ReferenceFE)
  @abstractmethod
end

"""
    Conformity(reffe::ReferenceFE) -> Conformity
"""
function Conformity(reffe::ReferenceFE)
  @abstractmethod
end

"""
    get_face_own_dofs(reffe::ReferenceFE,conf::Conformity) -> Vector{Vector{Int}}
"""
function get_face_own_dofs(reffe::ReferenceFE,conf::Conformity)
  @abstractmethod
end

function Conformity(reffe::ReferenceFE,conf::Conformity)
  conf
end

function Conformity(reffe::ReferenceFE,conf::Nothing)
  Conformity(reffe)
end

function Conformity(reffe::ReferenceFE,sym::Symbol)
  @abstractmethod
end

"""
    get_face_own_dofs(reffe::ReferenceFE) -> Vector{Vector{Int}}
"""
function get_face_own_dofs(reffe::ReferenceFE)
  conf = Conformity(reffe)
  get_face_own_dofs(reffe,conf)
end

function get_face_own_dofs(reffe::ReferenceFE,conf::Nothing)
  get_face_own_dofs(reffe)
end

function get_face_own_dofs(reffe::ReferenceFE,conf::L2Conformity)
  _get_face_own_dofs_l2(reffe)
end

function _get_face_own_dofs_l2(reffe::ReferenceFE)
  p = get_polytope(reffe)
  r = [Int[] for i in 1:num_faces(p)]
  r[end] = collect(1:num_dofs(reffe))
  r
end

"""
    get_face_own_dofs_permutations(reffe::ReferenceFE,conf::Conformity) -> Vector{Vector{Vector{Int}}}
"""
function get_face_own_dofs_permutations(reffe::ReferenceFE,conf::Conformity)
  face_own_dofs = get_face_own_dofs(reffe,conf)
  _trivial_face_own_dofs_permutations(face_own_dofs)
end

function _trivial_face_own_dofs_permutations(face_own_dofs)
  [ [collect(Int,1:length(dofs)),]  for dofs in face_own_dofs  ]
end

"""
    get_face_own_dofs_permutations(reffe::ReferenceFE) -> Vector{Vector{Vector{Int}}}
"""
function get_face_own_dofs_permutations(reffe::ReferenceFE)
  conf = Conformity(reffe)
  get_face_own_dofs_permutations(reffe,conf)
end

function get_face_own_dofs_permutations(reffe::ReferenceFE,conf::Nothing)
  get_face_own_dofs_permutations(reffe)
end

"""
    get_face_dofs(reffe::ReferenceFE) -> Vector{Vector{Int}}

Returns a vector of vector that, for each face, stores the
dofids in the closure of the face.
"""
function get_face_dofs(reffe::ReferenceFE)
  @abstractmethod
end

function get_dof_to_comp(reffe::ReferenceFE)
  fill(0,num_dofs(reffe))
end

# Push forward-related

abstract type PushForwardMap <: Map end

function evaluate!(cache,::PushForwardMap,v::AbstractVector{<:Field},phi::Field)
  @abstractmethod
end

PushForwardMap(::Type{<:ReferenceFE}) = IdentityPushForwardMap()

PushForwardMap(reffe::T) where T<:ReferenceFE = PushForwardMap(T)

struct IdentityPushForwardMap <: PushForwardMap end

function evaluate!(cache,::IdentityPushForwardMap,v::AbstractVector{<:Field},phi::Field)
  v
end

function lazy_map(::IdentityPushForwardMap,a::AbstractArray,b::AbstractArray)
  a
end

"""
"""
function get_shapefuns(reffe::ReferenceFE,phi::Field)
  PushForwardMap(reffe)(get_shapefuns(reffe),phi)
end

"""
"""
function get_dof_basis(reffe::ReferenceFE,phi::Field)
  get_dof_basis(reffe,phi,PushForwardMap(reffe))
end

function get_dof_basis(reffe::ReferenceFE,phi::Field,::IdentityPushForwardMap)
  get_dof_basis(reffe)
end

function get_dof_basis(reffe::ReferenceFE,phi::Field,::PushForwardMap)
  @abstractmethod
end

function lazy_map(::typeof(get_shapefuns),cell_reffe::AbstractArray,cell_map::AbstractArray)
  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
  ctype_ref_shapefuns = map(get_shapefuns,ctype_reffe)
  cell_ref_shapefuns = expand_cell_data(ctype_ref_shapefuns,cell_ctype)
  ctype_k = map(PushForwardMap,ctype_reffe)
  unique_ks = unique(ctype_k)
  if length(unique_ks) == 1
    k = first(unique_ks)
    lazy_map(k,cell_ref_shapefuns,cell_map)
  else
    T = return_type(get_shapefuns, testitem(cell_reffe), testitem(cell_map))
    lazy_map(get_shapefuns,T,cell_reffe,cell_map)
  end
end

function lazy_map(::typeof(get_dof_basis),cell_reffe::AbstractArray,cell_map::AbstractArray)
  ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
  if all( map(reffe->PushForwardMap(reffe)==IdentityPushForwardMap(),ctype_reffe) )
    ctype_dof_basis = map(get_dof_basis,ctype_reffe)
    expand_cell_data(ctype_dof_basis,cell_ctype)
  else
    T = return_type(get_dof_basis, testitem(cell_reffe), testitem(cell_map))
    lazy_map(get_dof_basis,T,cell_reffe,cell_map)
  end
end

"""
"""
function expand_cell_data(type_to_data, cell_to_type)
  CompressedArray(type_to_data,cell_to_type)
end

function expand_cell_data(type_to_data, cell_to_type::Fill)
  ncells = length(cell_to_type)
  @assert length(type_to_data) == 1 "Only one reference element expected"
  @assert cell_to_type.value == 1 "Only one type of reference element expected"
  data = first(type_to_data)
  Fill(data,ncells)
end

function compress_cell_data(cell_data::AbstractArray)
  @unreachable """\n
  The given cell data cannot be compressed. Describe your data with
  a CompressedArray or Fill array types.
  """
end

function compress_cell_data(a::CompressedArray)
  a.values, a.ptrs
end

function compress_cell_data(a::Fill)
  fill(a.value,1), Fill(1,length(a))
end

# Test

"""
    test_reference_fe(reffe::ReferenceFE{D}) where D

Test if the methods in the `ReferenceFE` interface are defined for the object `reffe`.
"""
function test_reference_fe(reffe::ReferenceFE{D}) where D
  conf = Conformity(reffe)
  @test isa(conf,Conformity)
  test_reference_fe(reffe,conf)
end

function test_reference_fe(reffe::ReferenceFE{D},conf::Conformity) where D
  @test D == num_dims(reffe)
  p = get_polytope(reffe)
  @test isa(p,Polytope{D})
  basis = get_prebasis(reffe)
  @test isa(basis,AbstractArray{<:Field})
  dofs = get_dof_basis(reffe)
  @test isa(dofs,AbstractArray{<:Dof})
  facedofs = get_face_own_dofs(reffe,conf)
  @test isa(facedofs,Vector{Vector{Int}})
  @test length(facedofs) == num_faces(p)
  facedofs_perms = get_face_own_dofs_permutations(reffe,conf)
  @test isa(facedofs_perms,Vector{Vector{Vector{Int}}})
  @test length(facedofs_perms) == num_faces(p)
  facedofs = get_face_dofs(reffe)
  @test isa(facedofs,Vector{Vector{Int}})
  @test length(facedofs) == num_faces(p)
  shapefuns = get_shapefuns(reffe)
  @test isa(shapefuns,AbstractVector{<:Field})
  ndofs = num_dofs(reffe)
  m = evaluate(dofs,basis)
  @test ndofs == size(m,1)
  @test ndofs == size(m,2)
end


"""
Constant of type `Int`  used to signal that a permutation is not valid.
"""
const INVALID_PERM = 0

# API

"""
    num_dims(::Type{<:ReferenceFE{D}}) where D
    num_dims(reffe::ReferenceFE{D}) where D

Returns `D`.
"""
num_dims(reffe::ReferenceFE) = num_dims(typeof(reffe))
num_dims(::Type{<:ReferenceFE{D}}) where D = D

"""
    num_cell_dims(::Type{<:ReferenceFE{D}}) where D
    num_cell_dims(reffe::ReferenceFE{D}) where D

Returns `D`.
"""
num_cell_dims(reffe::ReferenceFE) = num_dims(typeof(reffe))
num_cell_dims(::Type{<:ReferenceFE{D}}) where D = D

"""
    num_point_dims(::Type{<:ReferenceFE{D}}) where D
    num_point_dims(reffe::ReferenceFE{D}) where D

Returns `D`.
"""
num_point_dims(reffe::ReferenceFE) = num_dims(typeof(reffe))
num_point_dims(::Type{<:ReferenceFE{D}}) where D = D

"""
    num_faces(reffe::ReferenceFE)
    num_faces(reffe::ReferenceFE,d::Integer)
"""
num_faces(reffe::ReferenceFE) = num_faces(get_polytope(reffe))
num_faces(reffe::ReferenceFE,d::Integer) = num_faces(get_polytope(reffe),d)

"""
    num_vertices(reffe::ReferenceFE)
"""
num_vertices(reffe::ReferenceFE) = num_vertices(get_polytope(reffe))


"""
    num_edges(reffe::ReferenceFE)
"""
num_edges(reffe::ReferenceFE) = num_edges(get_polytope(reffe))


"""
    num_facets(reffe::ReferenceFE)
"""
num_facets(reffe::ReferenceFE) = num_facets(get_polytope(reffe))

"""
    get_face_own_dofs(reffe::ReferenceFE,conf::Conformity,d::Integer)
"""
function get_face_own_dofs(reffe::ReferenceFE,conf::Conformity,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_own_dofs(reffe,conf)[range]
end

"""
    get_face_own_dofs(reffe::ReferenceFE,d::Integer)
"""
function get_face_own_dofs(reffe::ReferenceFE,d::Integer)
  conf = Conformity(reffe)
  get_face_own_dofs(reffe,conf,d)
end

"""
    get_face_dofs(reffe::ReferenceFE,d::Integer)
"""
function get_face_dofs(reffe::ReferenceFE,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_dofs(reffe)[range]
end

"""
    get_face_own_dofs_permutations(reffe::ReferenceFE,conf::Conformity,d::Integer)
"""
function get_face_own_dofs_permutations(reffe::ReferenceFE,conf::Conformity,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_own_dofs_permutations(reffe,conf)[range]
end

"""
    get_face_own_dofs_permutations(reffe::ReferenceFE,d::Integer)
"""
function get_face_own_dofs_permutations(reffe::ReferenceFE,d::Integer)
  conf = Conformity(reffe)
  get_face_own_dofs_permutations(reffe,conf,d)
end

"""
    get_own_dofs_permutations(reffe::ReferenceFE,conf::Conformity)
"""
function get_own_dofs_permutations(reffe::ReferenceFE,conf::Conformity)
  n = num_faces(get_polytope(reffe))
  get_face_own_dofs_permutations(reffe,conf)[n]
end

"""
    get_own_dofs_permutations(reffe::ReferenceFE)
"""
function get_own_dofs_permutations(reffe::ReferenceFE)
  conf = Conformity(reffe)
  get_own_dofs_permutations(reffe,conf)
end

"""
    get_shapefuns(reffe::ReferenceFE) -> Field

Returns the basis of shape functions (i.e. the canonical basis)
associated with the reference FE. The result is encoded as a `Field` object.
"""
function get_shapefuns(reffe::ReferenceFE)
  dofs = get_dof_basis(reffe)
  prebasis = get_prebasis(reffe)
  compute_shapefuns(dofs,prebasis)
end

"""
    compute_shapefuns(dofs,prebasis)

Helper function used to compute the shape function basis
associated with the dof basis `dofs` and the basis `prebasis`.

It is equivalent to

    change = inv(evaluate(dofs,prebasis))
    linear_combination(change,prebasis) # i.e. transpose(change)*prebasis
"""
function compute_shapefuns(dofs,prebasis)
  change = inv(evaluate(dofs,prebasis))
  linear_combination(change,prebasis)
end

# Concrete implementation

"""
    struct GenericRefFE{T,D} <: ReferenceFE{D}
      # + private fields
    end

This type is a *materialization* of the `ReferenceFE` interface. That is, it is a
`struct` that stores the values of all abstract methods in the `ReferenceFE` interface.
This type is useful to build reference FEs from the underlying ingredients without
the need to create a new type.

Note that some fields in this `struct` are type unstable deliberately in order to simplify the
type signature. Don't access them in computationally expensive functions,
instead extract the required fields before and pass them to the computationally expensive function.
"""
struct GenericRefFE{T,D} <: ReferenceFE{D}
  ndofs::Int
  polytope::Polytope{D}
  prebasis::AbstractVector{<:Field}
  dofs::AbstractVector{<:Dof}
  conformity::Conformity
  metadata
  face_dofs::Vector{Vector{Int}}
  shapefuns::AbstractVector{<:Field}
  @doc """
        GenericRefFE{T}(
        ndofs::Int,
        polytope::Polytope{D},
        prebasis::AbstractVector{<:Field},
        dofs::AbstractVector{<:Dof},
        conformity::Conformity,
        metadata,
        face_dofs::Vector{Vector{Int}},
        shapefuns::AbstractVector{<:Field}=compute_shapefuns(dofs,prebasis)) where {T,D}

  Constructs a `GenericRefFE` object with the provided data.
  """
  function GenericRefFE{T}(
    ndofs::Int,
    polytope::Polytope{D},
    prebasis::AbstractVector{<:Field},
    dofs::AbstractVector{<:Dof},
    conformity::Conformity,
    metadata,
    face_dofs::Vector{Vector{Int}},
    shapefuns::AbstractVector{<:Field}=compute_shapefuns(dofs,prebasis)) where {T,D}

    new{T,D}(
      ndofs,
      polytope,
      prebasis,
      dofs,
      conformity,
      metadata,
      face_dofs,
      shapefuns)
  end
end

num_dofs(reffe::GenericRefFE) = reffe.ndofs

get_polytope(reffe::GenericRefFE) = reffe.polytope

get_prebasis(reffe::GenericRefFE) = reffe.prebasis

get_dof_basis(reffe::GenericRefFE) = reffe.dofs

Conformity(reffe::GenericRefFE) = reffe.conformity

get_face_dofs(reffe::GenericRefFE) = reffe.face_dofs

get_shapefuns(reffe::GenericRefFE) = reffe.shapefuns
