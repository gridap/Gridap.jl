
"""
    abstract type Conformity

`Conformity` subtypes are singletons representing a FE space conformity, that is
a Sobolev space that the global finite element space is a subset of.

The available conformities, along with the `Symbol`s that some high level
constructors accept, are:
- [`L2Conformity`](@ref); `:L2`
- [`GradConformity`](@ref) (alias `H1Conformity`); `:H1`, `:Hgrad`, `:C0`
- [`CurlConformity`](@ref); `:Hcurl`, `:HCurl`
- [`DivConformity`](@ref); `:Hdiv`, `:HDiv`
- [`CDConformity`](@ref)
"""
abstract type Conformity end

"""
    struct L2Conformity <: Conformity
"""
struct L2Conformity <: Conformity end

"""
    struct GradConformity <: Conformity
"""
struct GradConformity <: Conformity end
const H1Conformity = GradConformity

"""
    struct CurlConformity <: Conformity
"""
struct CurlConformity <: Conformity end

"""
    struct DivConformity <: Conformity
"""
struct DivConformity <: Conformity end


"""
    abstract type ReferenceFE{D} <: GridapType

Abstract type representing a Reference finite element. `D` is the underlying
coordinate space dimension. We follow the Ciarlet definition.
A reference finite element is defined by a polytope (cell topology), a basis of
an interpolation space of top of this polytope (denoted here as the prebasis),
and a basis of the dual of this space , i.e. the Degrees of Freedom (DoFs).
From this information one can compute the shape functions, i.e the canonical
basis w.r.t. the DoFs via `compute_shapefuns`. It is also possible to use
polynomial basis as shape functions, and compute the DoFs basis from a
prebasis of DoFs via `compute_dofs`. In the end, any `reffe::ReferenceFE` should
verify that

    evaluate(get_dof_basis(reffe), get_shapefuns(reffe))

is approximately the identity matrix.

In addition, this type encodes information about how the shape functions of the
reference finite element are "glued" with that of neighbors in order to build
globally conforming spaces, c.f. [`get_face_own_dofs`](@ref) and [`get_face_own_dofs_permutations`](@ref).

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

"""
    abstract type ReferenceFEName

Supertype for the reference finite element name singleton types, e.g. [`lagrangian`](@ref).
Instances are used to select a reference FE in the [`ReferenceFE`](@ref ReferenceFE(n::ReferenceFEName,a...;k...)) constructor.
"""
abstract type ReferenceFEName end

# Extensible factory function

"""
    ReferenceFE(name::ReferenceFEName[, T::Type], order;  kwargs...)
    ReferenceFE(name::ReferenceFEName[, T::Type], orders; kwargs...)
    ReferenceFE(F::Symbol, r, k, [, T::Type]; kwargs...)

Signatures defining a reference finite element (but yet unspecified cell polytope):
- an element [`name`](@ref ReferenceFEName), value type `T` and `order(s)`, or
- a FEEC family `F ∈ (:P⁻, :P, :Q⁻, :S)`, with polynomial order `r` and form order `k`.

# Arguments
- `T`: type of scalar components of the shape function values, `Float64` by default.
  For elements supporting Cartesian product space of a scalar one (e.g. [`lagrangian`](@ref)), this can be a tensor type like `VectorValue{2,Float64}`.
- `order::Int`: the polynomial order parameter, or
- `orders::NTuple{D,Int}`: a tuple of order per space dimension for anysotropic elements.

Keyword arguments are element specific, except
- `rotate_90::Bool=false`, set to true for div-conforming FEEC bases in 2D.


!!! warning
    This method only returns the tuple of its arguments, the actual Reference
    FE(s) is(are) only built once the polytope(s) is(are) known. See the other
    `ReferenceFE` methods or the FESpaces constructors.
"""
ReferenceFE(name::ReferenceFEName, args...; kwargs...) = (name, args, kwargs)
ReferenceFE(F::Symbol, args...; kwargs...) = (F, args, kwargs)

"""
    ReferenceFE(p::Polytope, args...; kwargs...)

Return the specified reference FE implemented on `p`.

The `args` and `kwargs` are the arguments of
[`ReferenceFE(::ReferenceFEName, ...; ...)`](@ref
ReferenceFE(::ReferenceFEName,a...;k...)) or [`ReferenceFE(F::Symbol, ...; ...)`](@ref
ReferenceFE(::Symbol,a...;k...)), first argument included.
"""
function ReferenceFE(p::Polytope, args...; kwargs...)
  @unreachable """\n
  Undefined factory function ReferenceFE for the given arguments:\n
  - args: $args,
  - kwargs: $kwargs.
  """
end
# Default component/value type
function ReferenceFE(p::Polytope, name::ReferenceFEName, order; kwargs...)
  ReferenceFE(p,name,Float64,order; kwargs...)
end
function ReferenceFE(p::Polytope,F::Symbol,r,k; kwargs...)
  ReferenceFE(p,F,r,k,Float64; kwargs...) # implemented in ExteriorCalculusRefFEs.jl
end

"""
    get_name(::ReferenceFE)
    get_name(::Type{ReferenceFE})

Returns the [ReferenceFEName](@ref) of the given reference FE.
"""
get_name(::Type{<:ReferenceFE})::ReferenceFEName = @abstractmethod
get_name(reffe::ReferenceFE) = get_name(typeof(reffe))

"""
    num_dofs(reffe::ReferenceFE) -> Int

Returns the number of DoFs, that is also the number of shape functions and the
dimension of the element's polynomial space.
"""
function num_dofs(reffe::ReferenceFE)
  @abstractmethod
end

"""
    get_polytope(reffe::ReferenceFE) -> Polytope

Returns the underlying [`Polytope`](@ref) object.
"""
function get_polytope(reffe::ReferenceFE)
  @abstractmethod
end

"""
    get_prebasis(reffe::ReferenceFE)

Returns the polynomial basis used as a pre-basis to compute the shape function,
which is a vector of `Field`s, often a [`PolynomialBasis`](@ref).
"""
function get_prebasis(reffe::ReferenceFE)
  @abstractmethod
end

"""
    get_order(reffe::ReferenceFE) -> Int

Returns the maximum polynomial order of shape function of `reffe`.
"""
get_order(reffe::ReferenceFE) = get_order(get_prebasis(reffe))

"""
    get_dof_basis(reffe::ReferenceFE)

Returns the underlying dof basis, which is a vector of `Dof`s, typically
[`LagrangianDofBasis{P,V}`](@ref) for FE using point evaluation based DoFs, or a
[`MomentBasedDofBasis`](@ref) using moment (face integral) based DoFs.
"""
function get_dof_basis(reffe::ReferenceFE)
  @abstractmethod
end

"""
    Conformity(reffe::ReferenceFE) -> Conformity

Return the default conformity of `reffe`.
"""
function Conformity(reffe::ReferenceFE)
  @abstractmethod
end

"""
    get_face_own_dofs(reffe::ReferenceFE[, conf::Conformity][, d::Int]) -> Vector{Vector{Int}}

Return a vector containing, for each face of `reffe`'s polytope, the indices of
the DoFs that belong to the interior of the face. This determines which DoFs
will be glued together in the global FE space: DoFs owned by a vertex or edge
shared by several physical polytopes will be glued/unified. Only the interior
owned DoFs (last face, i.e. the polytope itself) are not glued.

The ownership thus depends on the [`Conformity`](@ref) of the element: a `:L2`
conforming element is an element for which the polytope owns all DoFs such that
none are glued (e.g. for discontinuous Galerkin methods).

If `conf` is given, the ownership is computed for `conf` and not for `reffe`'s
conformity. If `d` is given, the returned vector only contains the data for the
`d`-dimensional faces of  `reffe`'s polytope.
"""
function get_face_own_dofs(reffe::ReferenceFE, conf::Conformity)
  @abstractmethod
end

function Conformity(reffe::ReferenceFE, conf::Conformity)
  conf
end

function Conformity(reffe::ReferenceFE, conf::Nothing)
  Conformity(reffe)
end

"""
    Conformity(reffe::ReferenceFE, conf::Symbol)

Return the [`Conformity`](@ref) that `conf` represents if `reffe` implements it,
or throws an `ErrorException` otherwise.

For example, if `reffe` is a Lagrangian refference FE with `H1Conformity`,
the function would return `L2Conformity()` and `H1Conformity()` for
respectively `conf=:L2` and `:H1` (because H¹ is in L²), but would error on
`conf=:Hcurl`.
"""
function Conformity(reffe::ReferenceFE, sym::Symbol)
  @abstractmethod
end

function get_face_own_dofs(reffe::ReferenceFE)
  conf = Conformity(reffe)
  get_face_own_dofs(reffe, conf)
end

function get_face_own_dofs(reffe::ReferenceFE, conf::Nothing)
  get_face_own_dofs(reffe)
end

function get_face_own_dofs(reffe::ReferenceFE, conf::L2Conformity)
  _get_face_own_dofs_l2(reffe)
end

function _get_face_own_dofs_l2(reffe::ReferenceFE)
  p = get_polytope(reffe)
  r = [Int[] for i in 1:num_faces(p)]
  r[end] = collect(1:num_dofs(reffe))
  r
end

"""
    get_face_own_dofs_permutations(
      reffe::ReferenceFE[, conf::Conformity][, d::Integer]) -> Vector{Vector{Vector{Int}}}

Return, for each `face` of `reffe`'s polytope, the permutation of the DoF indices
(that `face` owns) corresponding to each possible configuration of the `face`
once mapped in the physical space. These permutations correspond to the vertex
permutation returned by
[`get_face_vertex_permutations(get_polytope(reffe))`](@ref get_face_vertex_permutations).

`conf` can be passed to chose the conformity defining the DoFs ownership to the faces.
If `d` is given, the permutations are only returned for the `d`-dimensional faces.

In some cases, a valid permutation of the dofs for a corresponding
relabeling of the vertices does not exist (i.e., lagrangian FEs on `QUAD` with
anisotropic order). For these permutations the corresponding vector in
`get_face_own_dofs_permutations(reffe)` has to be filled with the constant value
[`INVALID_PERM`](@ref).
"""
function get_face_own_dofs_permutations(reffe::ReferenceFE, conf::Conformity)
  face_own_dofs = get_face_own_dofs(reffe, conf)
  _trivial_face_own_dofs_permutations(face_own_dofs)
end

function _trivial_face_own_dofs_permutations(face_own_dofs)
  [[collect(Int, 1:length(dofs)),] for dofs in face_own_dofs]
end

function get_face_own_dofs_permutations(reffe::ReferenceFE)
  conf = Conformity(reffe)
  get_face_own_dofs_permutations(reffe, conf)
end

function get_face_own_dofs_permutations(reffe::ReferenceFE, conf::Nothing)
  get_face_own_dofs_permutations(reffe)
end

"""
    get_face_dofs(reffe::ReferenceFE[, d::Integer]) -> Vector{Vector{Int}}

Returns a vector of vector that, for each face of `reffe`'s polytope, stores the
ids of the DoF in the closure of the face. The difference with
[`get_face_own_dofs`](@ref) is that this includes DoFs owned by the boundary
faces of each face.

If `d` is given, the returned vector only contains the data for the
`d`-dimensional faces of `reffe`'s polytope.
"""
function get_face_dofs(reffe::ReferenceFE)
  @abstractmethod
end

function get_dof_to_comp(reffe::ReferenceFE)
  fill(0, num_dofs(reffe))
end

"""
    expand_cell_data(id_to_data, cell_to_id) -> cell_to_data

Takes a (small) vector of possible datas `type_to_data`, and an array
`cell_to_type` of indices in `type_to_data`, and return a cell array
containing the data of each cell.

See also [`compress_cell_data`](@ref).
"""
function expand_cell_data(type_to_data, cell_to_type)
  CompressedArray(type_to_data, cell_to_type)
end

function expand_cell_data(type_to_data, cell_to_type::Fill)
  ncells = length(cell_to_type)
  @assert length(type_to_data) == 1 "Only one reference element expected"
  @assert cell_to_type.value == 1 "Only one type of reference element expected"
  data = first(type_to_data)
  Fill(data, ncells)
end

function expand_cell_data(type_to_data, cell_to_type::Base.OneTo)
  @assert length(type_to_data) == length(cell_to_type)
  type_to_data
end

function compress_cell_data(cell_data::AbstractArray)
  @unreachable """\n
  The given cell data cannot be compressed. Describe your data with
  a CompressedArray or Fill array types.
  """
end

"""
    compress_cell_data(cell_data::CompressedArray) -> id_to_data, cell_to_id
    compress_cell_data(cell_data::Fill) -> id_to_data, cell_to_id

Takes an array, iddentifies the unique values it contains, and returns the
vector of unique values and an array `cell_to_id` same shape as the input, but
containing the index of the data in `id_to_data` instead of the data.
See also [`CompressedArray`](@ref).

Use [`expand_cell_data`](@ref) to reverse the operation.
"""
function compress_cell_data(a::CompressedArray)
  a.values, a.ptrs
end

function compress_cell_data(a::Fill)
  fill(a.value, 1), Fill(1, length(a))
end

function compress_cell_data(a::Vector)
  a, Base.OneTo(length(a))
end

# Test

"""
    test_reference_fe(reffe::ReferenceFE{D}) where D

Test if the methods in the `ReferenceFE` interface are defined for the object `reffe`.
"""
function test_reference_fe(reffe::ReferenceFE{D}) where {D}
  conf = Conformity(reffe)
  @test isa(conf, Conformity)
  test_reference_fe(reffe, conf)
end

function test_reference_fe(reffe::ReferenceFE{D}, conf::Conformity) where {D}
  @test D == num_dims(reffe)
  p = get_polytope(reffe)
  @test isa(p, Polytope{D})
  basis = get_prebasis(reffe)
  @test isa(basis, AbstractArray{<:Field})
  dofs = get_dof_basis(reffe)
  @test isa(dofs, AbstractArray{<:Dof})
  facedofs = get_face_own_dofs(reffe, conf)
  @test isa(facedofs, Vector{Vector{Int}})
  @test length(facedofs) == num_faces(p)
  facedofs_perms = get_face_own_dofs_permutations(reffe, conf)
  @test isa(facedofs_perms, Vector{Vector{Vector{Int}}})
  @test length(facedofs_perms) == num_faces(p)
  facedofs = get_face_dofs(reffe)
  @test isa(facedofs, Vector{Vector{Int}})
  @test length(facedofs) == num_faces(p)
  shapefuns = get_shapefuns(reffe)
  @test isa(shapefuns, AbstractVector{<:Field})
  ndofs = num_dofs(reffe)
  m = evaluate(dofs, basis)
  @test ndofs == size(m, 1)
  @test ndofs == size(m, 2)
end


"""
Constant of type `Int`  used to signal that a permutation is not valid.
"""
const INVALID_PERM = 0

# API

"""
    num_dims(::Type{<:ReferenceFE{D}})
    num_dims(reffe::ReferenceFE{D})

Returns `D`.
"""
num_dims(reffe::ReferenceFE) = num_dims(typeof(reffe))
num_dims(::Type{<:ReferenceFE{D}}) where {D} = D

"""
    num_cell_dims(::Type{<:ReferenceFE{D}})
    num_cell_dims(reffe::ReferenceFE{D})

Returns `D`.
"""
num_cell_dims(reffe::ReferenceFE) = num_dims(typeof(reffe))
num_cell_dims(::Type{<:ReferenceFE{D}}) where {D} = D

"""
    num_point_dims(::Type{<:ReferenceFE{D}})
    num_point_dims(reffe::ReferenceFE{D})

Returns `D`.
"""
num_point_dims(reffe::ReferenceFE) = num_dims(typeof(reffe))
num_point_dims(::Type{<:ReferenceFE{D}}) where {D} = D

"""
    num_faces(reffe::ReferenceFE)
    num_faces(reffe::ReferenceFE, d::Integer)

Number of faces of `reffe`s polytope, counting the polytope itself.

If `d` is given, only the `d`-dimensional faces are counted.
"""
num_faces(reffe::ReferenceFE) = num_faces(get_polytope(reffe))
num_faces(reffe::ReferenceFE, d::Integer) = num_faces(get_polytope(reffe), d)

"""
    num_vertices(reffe::ReferenceFE)

Number of `0`-faces of `reffe`s polytope.
"""
num_vertices(reffe::ReferenceFE) = num_vertices(get_polytope(reffe))


"""
    num_edges(reffe::ReferenceFE)

Number of `1`-faces of `reffe`s polytope.
"""
num_edges(reffe::ReferenceFE) = num_edges(get_polytope(reffe))


"""
    num_facets(reffe::ReferenceFE)

Number of `(D-1)`-faces of `reffe`s polytope.
"""
num_facets(reffe::ReferenceFE) = num_facets(get_polytope(reffe))

function get_face_own_dofs(reffe::ReferenceFE, conf::Conformity, d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p, d)
  get_face_own_dofs(reffe, conf)[range]
end

function get_face_own_dofs(reffe::ReferenceFE, d::Integer)
  conf = Conformity(reffe)
  get_face_own_dofs(reffe, conf, d)
end

function get_face_dofs(reffe::ReferenceFE, d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p, d)
  get_face_dofs(reffe)[range]
end

function get_face_own_dofs_permutations(reffe::ReferenceFE, conf::Conformity, d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p, d)
  get_face_own_dofs_permutations(reffe, conf)[range]
end

function get_face_own_dofs_permutations(reffe::ReferenceFE, d::Integer)
  conf = Conformity(reffe)
  get_face_own_dofs_permutations(reffe, conf, d)
end

"""
    get_own_dofs_permutations(reffe::ReferenceFE)
    get_own_dofs_permutations(reffe::ReferenceFE, conf::Conformity)

Same as [`get_face_own_dofs_permutations`](@ref), but only return the last
permutations vector, that of DoFs owned by the cell (but not its boundary faces).
"""
function get_own_dofs_permutations(reffe::ReferenceFE, conf::Conformity)
  last(get_face_own_dofs_permutations(reffe, conf))
end

function get_own_dofs_permutations(reffe::ReferenceFE)
  conf = Conformity(reffe)
  get_own_dofs_permutations(reffe, conf)
end

"""
    get_shapefuns(reffe::ReferenceFE)

Returns the underlying shape function basis, which is an vector of `Field`s,
often a linear combination of a [`PolynomialBasis`](@ref) (the prebasis) and a
change of basis matrix.
"""
function get_shapefuns(reffe::ReferenceFE)
  dofs = get_dof_basis(reffe)
  prebasis = get_prebasis(reffe)
  compute_shapefuns(dofs, prebasis)
end

"""
    compute_shapefuns(dofs,prebasis)

Helper function used to compute the shape function basis
associated with the dof basis `dofs` and the basis `prebasis`.

It is equivalent to

    change = inv(evaluate(dofs,prebasis))
    linear_combination(change,prebasis) # i.e. transpose(change)*prebasis
"""
function compute_shapefuns(dofs, prebasis)
  change = inv(evaluate(dofs, prebasis))
  linear_combination(change, prebasis)
end

"""
    compute_dofs(predofs,shapefuns)

Helper function used to compute the dof basis
associated with the DoF prebasis `predofs` and the polynonial basis `shapefuns`.

It is equivalent to

    change = transpose(inv(evaluate(predofs,shapefuns)))
    linear_combination(change,predofs) # i.e. transpose(change)*predofs
"""
function compute_dofs(predofs, shapefuns)
  change = transpose(inv(evaluate(predofs, shapefuns)))
  linear_combination(change, predofs)
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
          shapefuns::AbstractVector{<:Field}=compute_shapefuns(dofs,prebasis)
        ) where {T,D}

  Constructor using a DoF basis and function space pre-basis.
  """
  function GenericRefFE{T}(
    ndofs::Int,
    polytope::Polytope{D},
    prebasis::AbstractVector{<:Field},
    dofs::AbstractVector{<:Dof},
    conformity::Conformity,
    metadata,
    face_dofs::Vector{Vector{Int}},
    shapefuns::AbstractVector{<:Field}=compute_shapefuns(dofs, prebasis)) where {T,D}

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
  @doc """
        GenericRefFE{T}(
          ndofs::Int,
          polytope::Polytope{D},
          prebasis::AbstractVector{<:Field},
          predofs::AbstractVector{<:Dof},
          conformity::Conformity,
          metadata,
          face_dofs::Vector{Vector{Int}},
          shapefuns::AbstractVector{<:Field},
          dofs::AbstractVector{<:Dof}=compute_dofs(predofs,shapefuns)
        ) where {T,D}

  Constructor using shape function basis and a DoF pre-basis.
  """
  function GenericRefFE{T}(
    ndofs::Int,
    polytope::Polytope{D},
    predofs::AbstractVector{<:Dof},
    conformity::Conformity,
    metadata,
    face_dofs::Vector{Vector{Int}},
    shapefuns::AbstractVector{<:Field},
    dofs::AbstractVector{<:Dof}=compute_dofs(predofs, shapefuns)) where {T,D}

    new{T,D}(
      ndofs,
      polytope,
      shapefuns,
      dofs,
      conformity,
      metadata,
      face_dofs,
      # Trick to be able to eval dofs af shapefuns in physical space
      # cf /test/FESpacesTests/PhysicalFESpacesTests.jl
      linear_combination(Eye{Int}(ndofs), shapefuns))
  end
end

get_name(::Type{<:GenericRefFE{Name}}) where Name = Name()

num_dofs(reffe::GenericRefFE) = reffe.ndofs

get_polytope(reffe::GenericRefFE) = reffe.polytope

get_prebasis(reffe::GenericRefFE) = reffe.prebasis

get_dof_basis(reffe::GenericRefFE) = reffe.dofs

Conformity(reffe::GenericRefFE) = reffe.conformity

get_face_dofs(reffe::GenericRefFE) = reffe.face_dofs

get_shapefuns(reffe::GenericRefFE) = reffe.shapefuns

get_metadata(reffe::GenericRefFE) = reffe.metadata

function ==(reffe1::GenericRefFE, reffe2::GenericRefFE)
  false
end

# TODO The hash is not consistent with this
function ==(reffe1::GenericRefFE{T,D}, reffe2::GenericRefFE{T,D}) where {T,D}
  t = true
  t = t && reffe1.ndofs       == reffe2.ndofs
  t = t && reffe1.polytope    == reffe2.polytope
  t = t && reffe1.prebasis    == reffe2.prebasis
  # Trick: only compare dofs OR shapefuns, one of those is a linear_combination
  # that does not implement ==
  t = t && ((reffe1.dofs      == reffe2.dofs)
        || (reffe1.shapefuns  == reffe2.shapefuns))
  t = t && reffe1.conformity  == reffe2.conformity
  t = t && reffe1.metadata    == reffe2.metadata
  t = t && reffe1.face_dofs   == reffe2.face_dofs
  t
end
