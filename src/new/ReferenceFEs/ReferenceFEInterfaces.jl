"""
    abstract type ReferenceFE{D}

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
- [`get_dofs(reffe::ReferenceFE)`](@ref)
- [`get_face_own_dofids(reffe::ReferenceFE)`](@ref)

and optionally these ones:

- [`ReferenceFE{N}(reffe::ReferenceFE,nfaceid::Integer) where N`](@ref)
- [`get_own_dofs_permutations(reffe::ReferenceFE)`](@ref)

"""
abstract type ReferenceFE{D} end

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

"""
    get_dofs(reffe::ReferenceFE) -> Dof

Returns the underlying dof basis encoded in a `Dof` object. 
"""
function get_dofs(reffe::ReferenceFE)
  @abstractmethod
end

"""
    get_face_own_dofids(reffe::ReferenceFE) -> Vector{Vector{Int}}
"""
function get_face_own_dofids(reffe::ReferenceFE)
  @abstractmethod
end

# optional

"""
    ReferenceFE{N}(reffe::ReferenceFE,nfaceid::Integer) where N

Returns a reference FE obtained by the restriction of the given one
to the face with `nfaceid` within dimension `N`.
"""
function ReferenceFE{N}(reffe::ReferenceFE,nfaceid::Integer) where N
  @abstractmethod
end

"""
    get_own_dofs_permutations(reffe::ReferenceFE) -> Vector{Vector{Int}}

Returns a vector of vectors indicating how the dofs owned by the reference
fe have to be relabeled when the vertices of the underlying polytope are relabeled
according the permutations in

    polytope = get_polytope(reffe)
    vertex_perms = get_vertex_permutations(polytope)

That is, if the vertices are relabeled with the permutation number `iperm`

    oldvertex_to_newvertex = vertex_perms[iperm]

The dofs are relabeled as

    olddof_to_newdof = get_own_dofs_permutations(reffe)[iperm]
 
Note that in some cases, a valid relabeling of the dofs
for a corresponding relabeling of the vertices does not exist (i.e., lagrangian FEs with
anisotropic order). For these permutations the corresponding vector in 
`get_own_dofs_permutations(reffe)` has to be filled with the constant value [`INVALID_PERM`](@ref).

# Examples

Relabeling of the dofs owned by the 4-th order segment

```jldoctest
using Gridap.ReferenceFEs

order = 4
seg5 = LagrangianRefFE(Float64,SEGMENT,order)

vertexperms = get_vertex_permutations(SEGMENT)
dofperms = get_own_dofs_permutations(seg5)

println(vertexperms)
println(dofperms)

# output
Array{Int64,1}[[1, 2], [2, 1]]
Array{Int64,1}[[1, 2, 3], [3, 2, 1]]

```

Note that when the segment is flipped (the second permutation) the first owned dof becomes the third 
owned dof, the second owned dof becomes the second one and the third becomes the first one as one would
expect.

Relabeling of the dofs owned by an anisotropic quadrilateral

```jldoctest
using Gridap.ReferenceFEs

orders = (2,3)
quad12 = LagrangianRefFE(Float64,QUAD,orders)

vertexperms = get_vertex_permutations(QUAD)
dofperms = get_own_dofs_permutations(quad12)

println(vertexperms)
println(dofperms)
println(INVALID_PERM)

# output
Array{Int64,1}[[1, 2, 3, 4], [1, 3, 2, 4], [2, 1, 4, 3], [2, 4, 1, 3], [3, 1, 4, 2], [3, 4, 1, 2], [4, 2, 3, 1], [4, 3, 2, 1]]
Array{Int64,1}[[1, 2], [0, 0], [1, 2], [0, 0], [0, 0], [2, 1], [0, 0], [2, 1]]
0
```

Note that not all permutations are valid in this case.

"""
function get_own_dofs_permutations(reffe::ReferenceFE)
  @abstractmethod
end

"""
Constant of type `Int`  used to signal that a permutation is not valid.
"""
const INVALID_PERM = 0


# API

"""
    get_shapefuns(reffe::ReferenceFE) -> Field

Returns the basis of shape functions (i.e. the canonical basis)
associated with the reference FE. The result is encoded as a `Field` object.
"""
function get_shapefuns(reffe::ReferenceFE)
  dofs = get_dofs(reffe)
  prebasis = get_prebasis(reffe)
  compute_shapefuns(dofs,prebasis)
end

"""
    compute_shapefuns(dofs,prebasis)

Helper function used to compute the shape function basis
associated with the dof basis `dofs` and the basis `prebasis`.

It is equivalent to

    change = inv(evaluate(dofs,prebasis))
    change_basis(prebasis,change)
"""
function compute_shapefuns(dofs,prebasis)
  change = inv(evaluate(dofs,prebasis))
  change_basis(prebasis,change)
end

num_dims(::Type{<:ReferenceFE{D}}) where D = D

"""
    num_dims(::Type{<:ReferenceFE{D}}) where D
    num_dims(reffe::ReferenceFE{D}) where D

Returns `D`.
"""
num_dims(reffe::ReferenceFE) = num_dims(typeof(reffe))

# Test

"""
    test_reference_fe(reffe::ReferenceFE{D};optional::Bool=false) where D

Test if the methods in the `ReferenceFE` interface are defined for the object `reffe`.
If `optional=false` (the default) only the mandatory methods are checked. Otherwise, the optional
methods are also tested.
"""
function test_reference_fe(reffe::ReferenceFE{D};optional::Bool=false) where D
  @test D == num_dims(reffe)
  p = get_polytope(reffe)
  @test isa(p,Polytope{D})
  basis = get_prebasis(reffe)
  @test isa(basis,Field)
  dofs = get_dofs(reffe)
  @test isa(dofs,Dof)
  facedofs = get_face_own_dofids(reffe)
  @test isa(facedofs,Vector{Vector{Int}})
  @test length(facedofs) == num_faces(p)
  shapefuns = get_shapefuns(reffe)
  @test isa(shapefuns,Field)
  ndofs = num_dofs(reffe)
  m = evaluate(dofs,basis)
  @test ndofs == size(m,1)
  @test ndofs == size(m,2)
  if optional
    dofperms = get_own_dofs_permutations(reffe)
    @test isa(dofperms,Vector{Vector{Int}})
    vertexperms = get_vertex_permutations(p)
    @test length(vertexperms) == length(dofperms)
    @test all( length.(dofperms) .== length(facedofs[end]) )
    for d in 0:D
      for iface in 1:num_faces(p,d)
        refface = ReferenceFE{d}(reffe,iface)
        @test isa(refface,ReferenceFE{d})
      end
    end
  end
end


# Concrete implementation

"""
    struct GenericRefFE{D} <: ReferenceFE{D}
      polytope::Polytope{D}
      prebasis::Field
      dofs::Dof
      face_own_dofids::Vector{Vector{Int}}
      shapefuns::Field
      ndofs::Int
      own_dofs_permutations::Vector{Vector{Int}}
      reffaces
    end

This type is a *materialization* of the `ReferenceFE` interface. That is, it is a 
`struct` that stores the values of all abstract methods in the `ReferenceFE` interface.
This type is useful to build reference FEs from the underlying ingredients without
the need to create a new type.

Note that some fields in this `struct` are type unstable deliberately in order to simplify the
type signature. Don't access them in computationally expensive functions,
instead extract the required fields before and pass them to the computationally expensive function.
"""
struct GenericRefFE{D} <: ReferenceFE{D}
    ndofs::Int
    polytope::Polytope{D}
    prebasis::Field
    dofs::Dof
    face_own_dofids::Vector{Vector{Int}}
    own_dofs_permutations::Vector{Vector{Int}}
    shapefuns::Field
    reffaces
  @doc """
      GenericRefFE(
        polytope::Polytope{D},
        prebasis::Field,
        dofs::Dof,
        face_own_dofids::Vector{Vector{Int}};
        shapefuns::Field,
        ndofs::Int,
        own_dofs_permutations::Vector{Vector{Int}},
        reffaces) where D

  Constructs a `GenericRefFE` object with the provided data. Positional arguments are mandatory.
  All Key-word ones are optional.
  """
  function GenericRefFE(
    polytope::Polytope{D},
    prebasis::Field,
    dofs::Dof,
    face_own_dofids::Vector{Vector{Int}};
    shapefuns::Field = compute_shapefuns(dofs,prebasis),
    ndofs::Int = size(evaluate(dofs,prebasis),1),
    own_dofs_permutations::Vector{Vector{Int}} = [ fill(INVALID_PERM,length(fi)) for fi in face_own_dofids],
    reffaces = nothing) where D

    new{D}(ndofs,polytope,prebasis,dofs,face_own_dofids,own_dofs_permutations,shapefuns,reffaces)
  end
end

num_dofs(reffe::GenericRefFE) = reffe.ndofs

get_polytope(reffe::GenericRefFE) = reffe.polytope

get_prebasis(reffe::GenericRefFE) = reffe.prebasis

get_dofs(reffe::GenericRefFE) = reffe.dofs

get_face_own_dofids(reffe::GenericRefFE) = reffe.face_own_dofids

get_own_dofs_permutations(reffe::GenericRefFE) = reffe.own_dofs_permutations

get_shapefuns(reffe::GenericRefFE) = reffe.shapefuns

function ReferenceFE{N}(reffe::GenericRefFE,iface::Integer) where N
  @assert reffe.reffaces != nothing "ReferenceFE cannot be provided. Make sure that you are using the keyword argument reffaces in the GenericRefFE constructor."
  reffe.reffaces[N+1][iface]
end

function ReferenceFE{D}(reffe::GenericRefFE{D},iface::Integer) where D
  @assert iface == 1 "Only one D-face"
  reffe
end

