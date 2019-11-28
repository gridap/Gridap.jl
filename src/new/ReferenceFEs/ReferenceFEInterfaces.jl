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
- [`get_dof_basis(reffe::ReferenceFE)`](@ref)
- [`get_face_own_dofs(reffe::ReferenceFE)`](@ref)

and optionally these ones:

- [`ReferenceFE{N}(reffe::ReferenceFE,nfaceid::Integer) where N`](@ref)
- [`get_own_dofs_permutations(reffe::ReferenceFE)`](@ref)
- [`(==)(a::ReferenceFE{D},b::ReferenceFE{D}) where D`](@ref)

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
    get_dof_basis(reffe::ReferenceFE) -> Dof

Returns the underlying dof basis encoded in a `Dof` object. 
"""
function get_dof_basis(reffe::ReferenceFE)
  @abstractmethod
end

"""
    get_face_own_dofs(reffe::ReferenceFE) -> Vector{Vector{Int}}
"""
function get_face_own_dofs(reffe::ReferenceFE)
  @abstractmethod
end

# optional

"""
    (==)(a::ReferenceFE{D},b::ReferenceFE{D}) where D

Returns `true` if the polytopes `a` and `b` are equivalent. Otherwise, it 
returns `false`.
"""
function (==)(a::ReferenceFE{D},b::ReferenceFE{D}) where D
  @abstractmethod
end

function (==)(a::ReferenceFE,b::ReferenceFE)
  false
end

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

"""
    get_face_dofs(reffe::ReferenceFE) -> Vector{Vector{Int}}

Returns a vector of vector that, for each face, stores the
dofids in the closure of the face.
"""
function get_face_dofs(reffe::ReferenceFE)
  polytope = get_polytope(reffe)
  D = num_dims(polytope)
  nfaces = num_faces(polytope)
  face_dofids = [Int[] for i in 1:nfaces]
  for d in 0:(D-1)
    _get_face_dofids_d!(face_dofids,Val{d}(),reffe,polytope)
  end
  face_dofids[end] = collect(1:num_dofs(reffe))
  face_dofids
end

function _get_face_dofids_d!(face_to_dofs,::Val{d},reffe,polytope) where d

  nface_to_own_dofs = get_face_own_dofs(reffe)
  nface_mface_to_nface = get_faces(polytope)
  offset = get_offset(polytope,d)

  for iface in 1:num_faces(polytope,d)
    face_reffe = ReferenceFE{d}(reffe,iface)
    face_polytope = get_polytope(face_reffe)
    mface_to_own_ldofs = get_face_own_dofs(face_reffe)
    dofs = zeros(Int,num_dofs(face_reffe))
    mface_to_nface = nface_mface_to_nface[iface+offset]
    for mface in 1:num_faces(face_polytope)
      nface = mface_to_nface[mface]
      own_dofs = nface_to_own_dofs[nface]
      own_ldofs = mface_to_own_ldofs[mface]
      dofs[own_ldofs] = own_dofs
    end
    face_to_dofs[iface+offset] = dofs
  end

end

"""
"""
function get_reffes(::Type{<:ReferenceFE{d}},reffe::ReferenceFE) where d
  ftype_to_reffe, _ = _compute_reffes_and_face_types(reffe,Val{d}())
  ftype_to_reffe
end

"""
"""
function get_face_types(::Type{<:ReferenceFE{d}},reffe::ReferenceFE) where d
  _, iface_to_ftype = _compute_reffes_and_face_types(reffe,Val{d}())
  iface_to_ftype
end

function _compute_reffes_and_face_types(reffe::ReferenceFE,::Val{d}) where d
  p = get_polytope(reffe)
  iface_to_reffe = [ ReferenceFE{d}(reffe,iface) for iface in 1:num_faces(p,d) ]
  _find_unique_with_indices(iface_to_reffe)
end

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
  dofs = get_dof_basis(reffe)
  @test isa(dofs,Dof)
  facedofs = get_face_own_dofs(reffe)
  @test isa(facedofs,Vector{Vector{Int}})
  @test length(facedofs) == num_faces(p)
  shapefuns = get_shapefuns(reffe)
  @test isa(shapefuns,Field)
  ndofs = num_dofs(reffe)
  m = evaluate(dofs,basis)
  @test ndofs == size(m,1)
  @test ndofs == size(m,2)
  if optional
    @test reffe == reffe
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
      face_own_dofs::Vector{Vector{Int}}
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
    face_own_dofs::Vector{Vector{Int}}
    own_dofs_permutations::Vector{Vector{Int}}
    shapefuns::Field
    reffaces
  @doc """
      GenericRefFE(
        polytope::Polytope{D},
        prebasis::Field,
        dofs::Dof,
        face_own_dofs::Vector{Vector{Int}};
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
    face_own_dofs::Vector{Vector{Int}};
    shapefuns::Field = compute_shapefuns(dofs,prebasis),
    ndofs::Int = size(evaluate(dofs,prebasis),1),
    own_dofs_permutations::Vector{Vector{Int}} = [ fill(INVALID_PERM,length(fi)) for fi in face_own_dofs],
    reffaces = nothing) where D

    new{D}(ndofs,polytope,prebasis,dofs,face_own_dofs,own_dofs_permutations,shapefuns,reffaces)
  end
end

num_dofs(reffe::GenericRefFE) = reffe.ndofs

get_polytope(reffe::GenericRefFE) = reffe.polytope

get_prebasis(reffe::GenericRefFE) = reffe.prebasis

get_dof_basis(reffe::GenericRefFE) = reffe.dofs

get_face_own_dofs(reffe::GenericRefFE) = reffe.face_own_dofs

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

"""
    abstract type NodalReferenceFE{D} <: ReferenceFE{D}

Abstract type representing a node-based reference FE. In addition
to all the information provided by a `ReferenceFE`, a nodal-based one
has information about the underlying nodes and relation between nodes
and DOFs.

The interface for this type is defined with the following methods

- [`get_node_coordinates(reffe::NodalReferenceFE)`](@ref)
- [`get_face_own_nodes(reffe::NodalReferenceFE)`](@ref)
- [`get_own_nodes_permutations(reffe::NodalReferenceFE)`](@ref)
- [`get_dof_to_node(reffe::NodalReferenceFE)`](@ref)
- [`get_dof_to_node(reffe::NodalReferenceFE)`](@ref)
- [`get_dof_to_comp(reffe::NodalReferenceFE)`](@ref)
- [`get_node_and_comp_to_dof(reffe::NodalReferenceFE)`](@ref)
"""
abstract type NodalReferenceFE{D} <: ReferenceFE{D} end

"""
    get_node_coordinates(reffe::NodalReferenceFE)
"""
function get_node_coordinates(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    get_face_own_nodes(reffe::NodalReferenceFE)
"""
function get_face_own_nodes(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    get_own_nodes_permutations(reffe::NodalReferenceFE)
"""
function get_own_nodes_permutations(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    get_dof_to_node(reffe::NodalReferenceFE)
"""
function get_dof_to_node(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    get_dof_to_comp(reffe::NodalReferenceFE)
"""
function get_dof_to_comp(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    get_node_and_comp_to_dof(reffe::NodalReferenceFE)
"""
function get_node_and_comp_to_dof(reffe::NodalReferenceFE)
  @abstractmethod
end

"""
    num_nodes(reffe::NodalReferenceFE)
"""
function num_nodes(reffe::NodalReferenceFE)
  length(get_node_coordinates(reffe))
end

"""
    get_face_nodes(reffe::NodalReferenceFE) -> Vector{Vector{Int}}

Returns a vector of vector that, for each face, stores the
nodeids in the closure of the face.
"""
function get_face_nodes(reffe::NodalReferenceFE)
  polytope = get_polytope(reffe)
  D = num_dims(polytope)
  nfaces = num_faces(polytope)
  face_nodeids = [Int[] for i in 1:nfaces]
  for d in 0:(D-1)
    _get_face_nodeids_d!(face_nodeids,Val{d}(),reffe,polytope)
  end
  face_nodeids[end] = collect(1:num_nodes(reffe))
  face_nodeids
end

function _get_face_nodeids_d!(face_to_nodes,::Val{d},reffe,polytope) where d

  nface_to_own_nodes = get_face_own_nodes(reffe)
  nface_mface_to_nface = get_faces(polytope)
  offset = get_offset(polytope,d)

  for iface in 1:num_faces(polytope,d)
    face_reffe = ReferenceFE{d}(reffe,iface)
    face_polytope = get_polytope(face_reffe)
    mface_to_own_lnodes = get_face_own_nodes(face_reffe)
    nodes = zeros(Int,num_nodes(face_reffe))
    mface_to_nface = nface_mface_to_nface[iface+offset]
    for mface in 1:num_faces(face_polytope)
      nface = mface_to_nface[mface]
      own_nodes = nface_to_own_nodes[nface]
      own_lnodes = mface_to_own_lnodes[mface]
      nodes[own_lnodes] = own_nodes
    end
    face_to_nodes[iface+offset] = nodes
  end

end

"""
    get_vertex_node(reffe::NodalReferenceFE) -> Vector{Int}
"""
function get_vertex_node(reffe::NodalReferenceFE)
  d = 0
  p = get_polytope(reffe)
  range = get_dimranges(p)[d+1]
  vertex_to_nodes = get_face_own_nodes(reffe)[range]
  map(first, vertex_to_nodes)
end

"""
    has_straight_faces(::NodalReferenceFE)

In the following sense:
vertices == nodes
"""
function has_straight_faces(reffe::NodalReferenceFE)
  p = get_polytope(reffe)
  r = true
  r = r && num_vertices(p) == num_nodes(reffe)
  r = r && get_vertex_node(reffe) == collect(1:num_nodes(reffe))
  r
end

"""
"""
function is_affine(reffe::NodalReferenceFE)
  p = get_polytope(reffe)
  has_straight_faces(reffe) && is_simplex(p)
end

"""
    test_nodal_reference_fe(reffe::NodalReferenceFE; optional::Bool=false)
"""
function test_nodal_reference_fe(reffe::NodalReferenceFE; optional::Bool=false)
  test_reference_fe(reffe,optional=optional)
  node_coordinates = get_node_coordinates(reffe)
  @test length(node_coordinates) == num_nodes(reffe)
  D = num_dims(reffe)
  @test isa(node_coordinates,Vector{<:Point{D}})
  face_own_nodes = get_face_own_nodes(reffe)
  @test isa(face_own_nodes,Vector{Vector{Int}})
  own_nodes_permutations = get_own_nodes_permutations(reffe)
  @test isa(own_nodes_permutations,Vector{Vector{Int}})
  dof_to_node = get_dof_to_node(reffe)
  @test isa(dof_to_node,Vector{Int})
  dof_to_comp = get_dof_to_comp(reffe)
  @test isa(dof_to_comp,Vector{Int})
  node_and_comp_to_dof = get_node_and_comp_to_dof(reffe)
  @test isa(node_and_comp_to_dof,Vector)
end

# IO

function Base.show(io::IO,p::ReferenceFE)
  print(io,nameof(typeof(p)))
end

