
# ---------------------------------------------------------------------------
# ConstraintLine
# ---------------------------------------------------------------------------

# One constrained (slave) DoF and its linear relation to master DoFs.
struct ConstraintLine{T<:Number}
  dof::Int
  masters::Vector{Int}
  coeffs::Vector{T}
  offset::T
end

function Base.show(io::IO, l::ConstraintLine)
  terms = join(["$(c)*x[$(d)]" for (d, c) in zip(l.masters, l.coeffs)], " + ")
  rhs = isempty(terms) ? "$(l.offset)" :
      iszero(l.offset) ? terms : "$terms + $(l.offset)"
  print(io, "x[$(l.dof)] = $rhs")
end

# ---------------------------------------------------------------------------
# ConstraintHandler
# ---------------------------------------------------------------------------

# DoF numbering: Like FESpaceWithLinearConstraints, we use unsigned DoF indices
# in ConstraintHandler. The function `add_constraint!` accepts signed DoF indices for
# user convenience, which get mapped into DOFs (unsigned).
# Then: 
# - `n_dofs`: Total number of DOFs in the space (free + dirichlet)
# - `n_fdofs`: Number of free DOFs (n_dofs >= n_fdofs)
mutable struct ConstraintHandler{T<:Number}
  n_dofs::Int
  n_fdofs::Int
  constraints::Vector{ConstraintLine{T}}
  dof_to_constraint::Vector{Int} # 0 if unconstrained
  closed::Bool
end

is_constrained(ch::ConstraintHandler, dof::Int) =
  !iszero(ch.dof_to_constraint[_dof_to_DOF(dof, ch.n_fdofs)])

num_constrained_dofs(ch::ConstraintHandler) = length(ch.constraints)
num_free_dofs(ch::ConstraintHandler)        = ch.n_dofs - num_constrained_dofs(ch)

function Base.show(io::IO, ch::ConstraintHandler{T}) where T
  status = ch.closed ? "closed" : "open"
  fdof_str = ch.n_fdofs < ch.n_dofs ? ", $(ch.n_fdofs) free" : ""
  print(io, "ConstraintHandler{$T}($(ch.n_dofs) dofs$fdof_str, " *
        "$(num_constrained_dofs(ch)) constrained, $status)")
end

"""
    ConstraintHandler(n_dofs, T=Float64)

Create an empty constraint handler for a space with `n_dofs` degrees of freedom.
All DOFs are treated as free (unsigned indices, no sign conversion).
`T` is the coefficient type.
"""
ConstraintHandler(n_dofs::Int, ::Type{T}=Float64) where T =
  ConstraintHandler{T}(n_dofs, n_dofs, ConstraintLine{T}[], zeros(Int, n_dofs), false)

"""
    ConstraintHandler(V::SingleFieldFESpace, T=Float64)

Create an empty constraint handler for the FE space `V`. DOF indices passed to
`add_constraint!` may use Gridap's signed convention (positive = free, negative =
Dirichlet); they are converted to unsigned indices internally.
`T` is the coefficient type.
"""
function ConstraintHandler(V::SingleFieldFESpace, ::Type{T}=Float64) where T
  n_fdofs = num_free_dofs(V)
  n_dofs  = n_fdofs + num_dirichlet_dofs(V)
  ConstraintHandler{T}(n_dofs, n_fdofs, ConstraintLine{T}[], zeros(Int, n_dofs), false)
end

"""
    ConstraintHandler(sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, V::FESpace, T=Float64)

Reconstruct an open `ConstraintHandler` from a Gridap Constructor-2 constraint
table triple defined on the `SingleFieldFESpace` `V`. This is the inverse of
`constraint_tables`.

The returned handler is open; call `close!` to chain-resolve and finalise it.
"""
function ConstraintHandler(
  sDOF_to_dof, sDOF_to_dofs::Table, sDOF_to_coeffs::Table,
  V::SingleFieldFESpace, ::Type{T}=Float64
) where T
  ch = ConstraintHandler(V, T)
  for (i, dof) in enumerate(sDOF_to_dof)
    add_constraint!(ch, dof, dataview(sDOF_to_dofs, i), dataview(sDOF_to_coeffs, i))
  end
  return ch
end

# ---------------------------------------------------------------------------
# Constraint addition
# ---------------------------------------------------------------------------

"""
    add_constraint!(ch, line; on_conflict=nothing)

Register a pre-built `ConstraintLine`. The other `add_constraint!` methods
delegate to this one.

When `on_conflict` is `nothing` (default), constraining an already-constrained
DoF raises an error. Otherwise `on_conflict(existing, incoming) → ConstraintLine`
is called and its return value replaces the stored line.
"""
function add_constraint!(
  ch::ConstraintHandler{T}, line::ConstraintLine{T};
  on_conflict = nothing
) where T
  @check !ch.closed "Cannot add constraints after close!()"
  @check 1 <= line.dof <= ch.n_dofs "DoF $(line.dof) out of range [1, $(ch.n_dofs)]"
  @check length(line.masters) == length(line.coeffs) "masters and coeffs must have the same length"
  for master in line.masters
    @check 1 <= master <= ch.n_dofs "Master DoF $master out of range [1, $(ch.n_dofs)]"
    @check master != line.dof "Self-referential constraint on DoF $(line.dof)"
  end

  dof = line.dof
  if !is_constrained(ch, dof)
    push!(ch.constraints, line)
    ch.dof_to_constraint[dof] = length(ch.constraints)
  elseif !isnothing(on_conflict)
    idx = ch.dof_to_constraint[dof]
    ch.constraints[idx] = on_conflict(ch.constraints[idx], line)
  else
    @notimplemented """ 
      Constraint conflict on DoF $dof:
        - existing: $(ch.constraints[ch.dof_to_constraint[dof]])
        - incoming: $line
      To resolve conflicts, provide an `on_conflict` function to `add_constraint!`
    """
  end
  return ch
end

"""
    add_constraint!(ch, dof, masters, coeffs, offset=zero(T); on_conflict=nothing)

Register the linear constraint  x[dof] = Σ coeffs[i] * x[masters[i]] + offset.

Both `masters` and `coeffs` are copied internally, so the caller can safely reuse
or modify them after the call.
"""
function add_constraint!(
  ch::ConstraintHandler{T}, dof::Int,
  masters::AbstractVector{<:Integer},
  coeffs::AbstractVector{<:Number},
  offset::Number = zero(T);
  kwargs...
) where T
  @check !ch.closed "Cannot add constraints to a closed handler"
  @check length(masters) == length(coeffs)
  DOF = _dof_to_DOF(dof, ch.n_fdofs)
  _masters = _dof_to_DOF!(collect(Int, masters), ch.n_fdofs)
  line = ConstraintLine{T}(DOF, _masters, collect(T, coeffs), T(offset))
  add_constraint!(ch, line; kwargs...)
  return ch
end

function add_constraint!(
  ch::ConstraintHandler{T},
  dofs::AbstractVector{<:Integer},
  masters::AbstractVector{<:Integer},
  coeffs::AbstractMatrix{<:Number},
  offsets::AbstractVector{<:Number} = zeros(T, size(coeffs, 1));
  kwargs...
) where T
  n_s, n_m = length(dofs), length(masters)
  @check size(coeffs, 1) == n_s "coeffs has $(size(coeffs,1)) rows but dofs has $n_s entries"
  @check size(coeffs, 2) == n_m "coeffs has $(size(coeffs,2)) cols but masters has $(length(masters)) entries"
  @check length(offsets) == n_s "offsets has $(length(offsets)) entries but dofs has $n_s entries"
  @inbounds for j in 1:n_s
    add_constraint!(
      ch, dofs[j], masters, view(coeffs, j, :), offsets[j]; kwargs...
    )
  end
  return ch
end

"""
    add_constraint!(ch, dof, offset)

Register a Dirichlet constraint  x[dof] = offset  (no master DoFs).
"""
function add_constraint!(ch::ConstraintHandler{T}, dof::Int, offset::Number; kwargs...) where T
  DOF = _dof_to_DOF(dof, ch.n_fdofs)
  add_constraint!(ch, ConstraintLine{T}(DOF, Int[], T[], T(offset)); kwargs...)
end

"""
    set_offset!(ch, dof, value)

Update the offset of an already-registered constraint on `dof`.
"""
function set_offset!(ch::ConstraintHandler{T}, dof::Int, value::Number) where T
  @check !ch.closed "Cannot modify constraints after close!()"
  @check is_constrained(ch, dof) "DoF $dof has no constraint to update"
  idx = ch.dof_to_constraint[_dof_to_DOF(dof, ch.n_fdofs)]
  line = ch.constraints[idx]
  ch.constraints[idx] = ConstraintLine{T}(line.dof, line.masters, line.coeffs, T(value))
  return ch
end

"""
    merge!(a, b; on_conflict=nothing)

Absorb all constraints from `b` into `a`. Both handlers must be open.
Conflict resolution is optional and controlled by the user-defined function `on_conflict`,
which is called when the same DoF is constrained in both handlers.
Its signature should be

  `on_conflict(existing::ConstraintLine, incoming::ConstraintLine) → ConstraintLine`

Example resolvers:
  - keep existing:  `(e, _) -> e`
  - keep incoming:  `(_, i) -> i`
  - average:        `(e, i) -> ConstraintLine(e.dof, vcat(e.masters, i.masters), vcat(e.coeffs, i.coeffs) .* 0.5, 0.0)`

Chain resolution and cycle detection are deferred to `close!(a)`.
"""
function Base.merge!(a::ConstraintHandler{T}, b::ConstraintHandler{T}; on_conflict=nothing) where T
  @check !a.closed && !b.closed "Both handlers must be open to merge"
  @check a.n_dofs == b.n_dofs "Handlers have incompatible sizes: $(a.n_dofs) vs $(b.n_dofs)"
  for line in b.constraints
    add_constraint!(a, line; on_conflict)
  end
  return a
end

"""
    close!(ch; tol=1e-13)

Finalize the constraint handler:
1. Sort constraints by DoF index; rebuild dof_to_constraint.
2. Compute a topological ordering of the slave→slave dependency graph,
   erroring on cycles.
3. Resolve constraint chains by single-pass substitution in topological order.
4. Merge duplicate master entries, strip near-zero coefficients.
"""
function close!(ch::ConstraintHandler{T}; tol::Float64=1e-13) where T
  ch.closed && return ch

  # Step 1: sort, rebuild dof_to_constraint.
  sort!(ch.constraints, by=l->l.dof)
  fill!(ch.dof_to_constraint, 0)
  for (i, line) in enumerate(ch.constraints)
    ch.dof_to_constraint[line.dof] = i
  end

  # Step 2+3: topological ordering, then single-pass chain resolution.
  # By the time slave `s` is processed, every slave `t` it depends on has
  # already been processed and fully resolved (no slave masters remain).
  order, _, has_chains = compute_topological_ordering(ch)
  dof_to_sDOF = ch.dof_to_constraint
  if has_chains
    masters, coeffs = Int[], T[] # reusable buffers
    for s in order
      line = ch.constraints[s]
      any(m -> dof_to_sDOF[m] != 0, line.masters) || continue

      empty!(masters); empty!(coeffs)
      new_offset = line.offset
      for (master, coeff) in zip(line.masters, line.coeffs)
        t = dof_to_sDOF[master]
        if t != 0
          sub = ch.constraints[t]
          append!(masters, sub.masters)
          append!(coeffs, coeff .* sub.coeffs)
          new_offset += coeff * sub.offset
        else
          push!(masters, master); push!(coeffs, coeff)
        end
      end
      resize!(line.masters, length(masters)); copyto!(line.masters, masters)
      resize!(line.coeffs,  length(coeffs)); copyto!(line.coeffs,  coeffs)
      ch.constraints[s] = ConstraintLine{T}(line.dof, line.masters, line.coeffs, new_offset)
    end
  end

  # Step 4: merge duplicate masters, strip near-zeros.
  # Masters and coeffs are modified in-place — no new ConstraintLine needed.
  z = zero(T)
  m_to_c = OrderedDict{Int,T}()
  for line in ch.constraints
    empty!(m_to_c)
    for (master, coeff) in zip(line.masters, line.coeffs)
      abs(coeff) < tol && continue
      m_to_c[master] = get(m_to_c, master, z) + coeff
    end
    n_new = length(m_to_c)
    resize!(line.masters, n_new)
    resize!(line.coeffs,  n_new)
    for (i, (master, coeff)) in enumerate(m_to_c)
      line.masters[i] = master
      line.coeffs[i] = coeff
    end
  end

  ch.closed = true
  return ch
end

"""
    compute_constraint_dag(ch::ConstraintHandler)
        → (dag::Table{Int}, in_degree::Vector{Int})

Build the slave→slave dependency graph of `ch`: slave `s` depends on slave
`t` if `t`'s DoF appears as a master of `s`. `dag[t]` lists `t`'s dependents
(CSR `Table`), for use by [`compute_topological_ordering`](@ref).
See the core method (in `FESpacesWithLinearConstraints.jl`) for details.

Requires `ch.constraints` to already be sorted by DoF (i.e. `ch.dof_to_constraint`
up to date), as is the case after step 1 of `close!`.
"""
function compute_constraint_dag(ch::ConstraintHandler)
  n_slaves    = length(ch.constraints)
  dof_to_sDOF = ch.dof_to_constraint

  # Build dependency graph (CSR): edge t → s means slave s depends on slave t.
  dep_ptrs  = zeros(Int32, n_slaves + 1)
  in_degree = zeros(Int, n_slaves)
  for s in 1:n_slaves
    for master in ch.constraints[s].masters
      t = dof_to_sDOF[master]
      if t != 0
        dep_ptrs[t+1] += Int32(1)
        in_degree[s]  += 1
      end
    end
  end

  Arrays.length_to_ptrs!(dep_ptrs)
  dep_data = Vector{Int}(undef, dep_ptrs[end] - 1)
  fill_pos = copy(dep_ptrs)
  for s in 1:n_slaves
    for master in ch.constraints[s].masters
      t = dof_to_sDOF[master]
      if t != 0
        dep_data[fill_pos[t]] = s
        fill_pos[t] += 1
      end
    end
  end

  dag = Table(dep_data, dep_ptrs)
  return dag, in_degree
end

"""
    get_matrix(ch) → (C, g)

Return the sparse prolongation matrix `C` (n × n_free) and offset vector `g`
(length n) such that the full solution vector satisfies `x = C * u + g`,
where `u` is the reduced (free-dof) solution.

Columns of `C` correspond to free DoFs in natural index order.
Free-DoF rows are identity rows; slave-DoF rows carry the constraint coefficients.
"""
function Algebra.get_matrix(ch::ConstraintHandler{T}) where T
  @check ch.closed "The ConstraintHandler must be closed before building the constraint matrix."
  n      = ch.n_dofs
  n_free = num_free_dofs(ch)

  # Map each free DoF to its column index in C.
  col = 0
  n_entries = 0
  dof_to_col = zeros(Int, n)
  for i in 1:n
    if !is_constrained(ch, i)
      col += 1
      dof_to_col[i] = col
      n_entries += 1
    else
      line = ch.constraints[ch.dof_to_constraint[i]]
      n_entries += length(line.masters)
    end
  end

  I = Vector{Int}(undef,n_entries)
  J = Vector{Int}(undef,n_entries)
  V = Vector{T}(undef,n_entries)
  g = zeros(T, n)

  k = 1
  for i in 1:n
    if !is_constrained(ch, i)
      I[k] = i
      J[k] = dof_to_col[i]
      V[k] = one(T)
      k += 1
    else
      line = ch.constraints[ch.dof_to_constraint[i]]
      g[i] = line.offset
      for (master, coeff) in zip(line.masters, line.coeffs)
        I[k] = i
        J[k] = dof_to_col[master]
        V[k] = coeff
        k += 1
      end
    end
  end

  C = sparse(I, J, V, n, n_free)
  return C, g
end

"""
    apply_constraints(A, b, ch) → (Ã, b̃)

Return the reduced system `Ã = C'AC` and `b̃ = C'(b - Ag)` obtained by
eliminating all constrained DoFs via the prolongation matrix `C` and offset
vector `g` from `constraint_matrix(ch)`.

Works with any matrix type; sparsity is preserved when `A` is sparse.
After solving `Ã * u = b̃`, recover the full solution with `x = C * u + g`.
"""
function apply_constraints(A, b::AbstractVector, ch::ConstraintHandler)
  @check ch.closed "Must call close!(ch) before condense"
  C, g = get_matrix(ch)
  Ã = C' * A * C
  b̃ = C' * (b - A * g)
  return Ã, b̃
end

"""
    apply_constraints!(x, ch)

Overwrite `x[dof]` for every constrained DoF using the stored constraints:

  x[dof] = offset + Σ_j coefficient_j * x[master_j]

After `close!()`, all masters are free, so the iteration order is immaterial.
"""
function apply_constraints!(x::AbstractVector, ch::ConstraintHandler)
  @check ch.closed "Must call close!() before apply_constraints!()"
  @check length(x) == ch.n_dofs "Vector length $(length(x)) does not match n_dofs $(ch.n_dofs)"
  for line in ch.constraints
    val = line.offset
    for (master, coeff) in zip(line.masters, line.coeffs)
      val += coeff * x[master]
    end
    x[line.dof] = val
  end
  return x
end

# ---------------------------------------------------------------------------
# FESpaceWithLinearConstraints from a ConstraintHandler
# ---------------------------------------------------------------------------

function constraint_tables(V::SingleFieldFESpace, ch::ConstraintHandler{T}) where T
  @check ch.closed "Must call close!(ch) before building constraint tables"

  n_fdofs = num_free_dofs(V)
  n_s     = length(ch.constraints)

  ptrs = Vector{Int}(undef, n_s + 1)
  ptrs[1] = 1
  for (i, line) in enumerate(ch.constraints)
    ptrs[i+1] = ptrs[i] + length(line.masters)
  end
  nnz          = ptrs[end] - 1
  master_dofs  = Vector{Int}(undef, nnz)
  coefficients = Vector{T}(undef, nnz)
  for (i, line) in enumerate(ch.constraints)
    for j in eachindex(line.masters)
      master_dofs[ ptrs[i] + j - 1] = _DOF_to_dof(line.masters[j], n_fdofs)
      coefficients[ptrs[i] + j - 1] = line.coeffs[j]
    end
  end

  sDOF_to_dof     = [_DOF_to_dof(line.dof, n_fdofs) for line in ch.constraints]
  sDOF_to_dofs    = Table(master_dofs, ptrs)
  sDOF_to_coeffs  = Table(coefficients, copy(ptrs))
  sDOF_to_offsets = T[line.offset for line in ch.constraints]

  return sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, sDOF_to_offsets
end

"""
    FESpaceWithLinearConstraints(V::SingleFieldFESpace, ch::ConstraintHandler)

Build a `FESpaceWithLinearConstraints` from `V` and a closed `ConstraintHandler`.
"""
function FESpaceWithLinearConstraints(V::SingleFieldFESpace, ch::ConstraintHandler)
  sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, sDOF_to_offsets = constraint_tables(V, ch)
  return FESpaceWithLinearConstraints(sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, V; sDOF_to_offsets)
end

# ---------------------------------------------------------------------------
# ConstraintHandler — cell-wise constraint arrays
# ---------------------------------------------------------------------------

"""
    ConstraintHandler(n_dofs, i_to_slave_ids, i_to_master_ids, i_to_coeffs, i_to_offsets; n_fdofs, on_conflict)

Build a `ConstraintHandler` from index-wise constraint data. Each index `i` contributes
one or more affine constraints. Two element shapes are supported:

**Multiple constraints per index** (`i_to_coeffs[i]` is a matrix `n_s × n_m`):

    x[i_to_slave_ids[i][j]] = Σ_k i_to_coeffs[i][j,k] * x[i_to_master_ids[i][k]] + i_to_offsets[i][j]

**Single constraint per index** (`i_to_coeffs[i]` is a vector `n_m`):

    x[i_to_slave_ids[i]] = i_to_coeffs[i] ⋅ x[i_to_master_ids[i]] + i_to_offsets[i]

Arrays may be lazy or concrete; iterated via the Gridap `array_cache`/`getindex!` protocol.
`T` is deduced as `promote_type(eltype(eltype(i_to_coeffs)), eltype(eltype(i_to_offsets)))`.

# Arguments
- `n_dofs`: total number of DOFs in the handler
- `i_to_slave_ids`: slave DOF index (scalar) or vector of indices, one per entry
- `i_to_master_ids`: master DOF index vector, one per entry
- `i_to_coeffs`: coefficient matrix or vector, one per entry
- `i_to_offsets`: offset vector or scalar, one per entry
- `n_fdofs`: free DOFs for signed→unsigned conversion (default: `n_dofs`, no conversion)
- `on_conflict`: called as `on_conflict(existing, incoming) → ConstraintLine`
  when the same DOF is constrained twice; default `nothing` raises an error
"""
function ConstraintHandler(
  n_dofs::Integer,
  i_to_slave_ids,
  i_to_master_ids,
  i_to_coeffs,
  i_to_offsets;
  n_fdofs::Int = Int(n_dofs),
  on_conflict = nothing
)
  T = promote_type(
    eltype(eltype(i_to_coeffs)),
    eltype(eltype(i_to_offsets))
  )
  ch = ConstraintHandler{T}(Int(n_dofs), n_fdofs, ConstraintLine{T}[], zeros(Int, Int(n_dofs)), false)

  cs = array_cache(i_to_slave_ids)
  cm = array_cache(i_to_master_ids)
  cc = array_cache(i_to_coeffs)
  co = array_cache(i_to_offsets)

  for i in eachindex(i_to_slave_ids)
    slaves  = getindex!(cs, i_to_slave_ids,  i)
    masters = getindex!(cm, i_to_master_ids, i)
    coeffs  = getindex!(cc, i_to_coeffs,     i)
    offsets = getindex!(co, i_to_offsets,    i)
    add_constraint!(ch, slaves, masters, coeffs, offsets; on_conflict)
  end

  return ch
end

# ---------------------------------------------------------------------------
# Automatic slave selection
# ---------------------------------------------------------------------------

"""
    select_constraint_slaves(n_dofs, cell_dof_ids, cell_coeffs, cell_offsets[, cell_nc[, T]]; atol)
    → (slave_ids, masters::Table, coeffs_table::Table, offsets_vec)

Given cell-wise constraint data, automatically selects one slave DOF per constraint
and returns the resulting tables.

All cells must use the same element shape — either uniformly `(matrix, vector)` or
uniformly `(vector, scalar)`:

- **Multiple constraints per cell**: `cell_coeffs[i]` is a matrix `n_s × n_m`,
  `cell_offsets[i]` is a vector of length `n_s`.
- **Single constraint per cell**: `cell_coeffs[i]` is a vector of length `n_m`,
  `cell_offsets[i]` is a scalar.

Supports both unsigned (all-positive) and Gridap's signed (free: positive, Dirichlet:
negative) DOF conventions. Pass `n_fdofs` to enable signed input.

**Slave selection.** For each constraint, the slave is the DOF with the largest
`|coefficient|`, breaking ties by choosing the DOF appearing in the fewest
remaining constraints (lower connectivity).

# Arguments
- `n_dofs`: total number of DOFs (upper bound on any DOF index in `cell_dof_ids`)
- `cell_dof_ids`: array of DOF ID vectors, one per cell (lazy or concrete)
- `cell_coeffs`: array of coefficient matrices or vectors, one per cell
- `cell_offsets`: array of offset vectors or scalars, one per cell
- `cell_nc`: constraints per cell (default: inferred from `cell_offsets`)
- `T`: numeric type (default `Float64`)
- `n_fdofs`: number of free DOFs; when `n_fdofs < n_dofs`, signed Gridap DOF IDs are
  accepted (default: `n_dofs`)
- `atol`: coefficients with `|a| ≤ atol` are treated as zero (default `eps(T)`)

# Returns
- `slave_ids::Vector{Int}`: selected slave DOF for each constraint
- `masters::Table{Int}`: master DOF IDs per constraint (CSR)
- `coeff_table::Table{T}`: normalised coefficients `−aᵢ/aₛ` per constraint (CSR)
- `offsets_vec::Vector{T}`: normalised affine offsets `bᵢ/aₛ` per constraint
"""
function select_constraint_slaves(
  n_dofs::Int,
  cell_dof_ids, cell_coeffs, cell_offsets,
  cell_nc = Arrays.lazy_collect(lazy_map(Base.Fix2(size, 1), cell_offsets)),
  ::Type{T} = Float64;
  n_fdofs::Int = n_dofs,
  atol::Real = eps(T)
) where T

  n_cells = length(cell_dof_ids)

  cd  = array_cache(cell_dof_ids)
  cc  = array_cache(cell_coeffs)
  co  = array_cache(cell_offsets)

  # Preliminary work: count total constraints and connectivity
  n_const = 0
  n_data  = 0
  dof_to_nconstraints = zeros(Int, n_dofs)
  for i in 1:n_cells
    dofs = getindex!(cd, cell_dof_ids, i)
    nci  = cell_nc[i]
    for d in dofs; dof_to_nconstraints[_dof_to_DOF(d, n_fdofs)] += nci; end
    n_const += nci
    n_data  += length(dofs) * nci
  end

  ptrs = zeros(Int32, n_const + 1)
  slaves       = Vector{Int}(undef, n_const)
  masters_data = Vector{Int}(undef, n_data)
  coeffs_data  = Vector{T}(undef, n_data)
  offsets      = Vector{T}(undef, n_const)
  is_slave     = zeros(Bool, n_dofs)

  k = 0
  p = 1
  for cell in 1:n_cells

    dofs = getindex!(cd, cell_dof_ids, cell)
    A    = getindex!(cc, cell_coeffs, cell)
    b    = getindex!(co, cell_offsets, cell)
    nd   = length(dofs)
    nci  = cell_nc[cell]

    IA = LinearIndices((nci, nd))
    Ib = LinearIndices((nci,))

    for i in 1:nci
      k += 1

      # Find slave: max |coeff|, break ties by lower connectivity
      best_j = 0; best_a = zero(T); best_conn = typemax(Int)
      for j in 1:nd
        DOF  = _dof_to_DOF(dofs[j], n_fdofs)
        a    = T(A[IA[i, j]])
        conn = dof_to_nconstraints[DOF]
        if abs(a) > atol && !is_slave[DOF]
          if abs(a) > abs(best_a) || (abs(a) == abs(best_a) && conn < best_conn)
            best_j = j; best_a = a; best_conn = conn
          end
        end
      end
      @check best_j > 0 "No slave found for constraint $k (cell $cell, row $i)"

      slaves[k]  = dofs[best_j]          # preserve input convention
      offsets[k] = T(b[Ib[i]]) / best_a
      is_slave[_dof_to_DOF(dofs[best_j], n_fdofs)] = true

      cnt = 0
      for j in 1:nd
        DOF = _dof_to_DOF(dofs[j], n_fdofs)
        a   = T(A[IA[i, j]])
        if j != best_j && abs(a) > atol
          masters_data[p] = dofs[j]; coeffs_data[p] = -a / best_a  # preserve input convention
          p += 1; cnt += 1
        end
        dof_to_nconstraints[DOF] -= 1
      end
      ptrs[k+1] = Int32(cnt)
    end
  end

  length_to_ptrs!(ptrs)
  resize!(masters_data, p - 1)
  resize!(coeffs_data,  p - 1)

  return slaves, Table(masters_data, ptrs), Table(coeffs_data, ptrs), offsets
end

# ---------------------------------------------------------------------------
# Internal helper — shared tail for all four weak-form constructors.
# Applies local_map if provided, runs slave selection, and builds the handler.
# ---------------------------------------------------------------------------

function _build_from_cells(
  space::FESpace, dof_ids, cell_coeffs, cell_offsets, cell_nc, ::Type{T};
  local_map  = nothing,
  on_conflict = nothing
) where T
  if !isnothing(local_map)
    cell_coeffs, cell_offsets = unpair_arrays(lazy_map(local_map, cell_coeffs, cell_offsets))
  end
  n_fdofs = num_free_dofs(space)
  n_dofs  = n_fdofs + num_dirichlet_dofs(space)
  slave_ids, masters, coeffs_t, offsets_v = select_constraint_slaves(
    n_dofs, dof_ids, cell_coeffs, cell_offsets, cell_nc, T; n_fdofs
  )
  return ConstraintHandler(n_dofs, slave_ids, masters, coeffs_t, offsets_v; n_fdofs, on_conflict)
end

# ---------------------------------------------------------------------------
# Constructor 1 — cell, bilinear: lhs(u,v), rhs(v)
# ---------------------------------------------------------------------------

"""
    ConstraintHandler(space, space_test, lhs, rhs; local_map, on_conflict, T)

Build a `ConstraintHandler` from a bilinear form `lhs(u,v)` and linear form `rhs(v)`.
One slave DOF is selected per constraint (= per test DOF per cell) via
`select_constraint_slaves`.

# Arguments
- `space`: trial FE space (master candidates)
- `space_test`: test FE space (one constraint per test DOF per cell)
- `lhs`: `(u,v) → DomainContribution` — cell coefficient matrices
- `rhs`: `v → DomainContribution` — cell offset vectors
- `local_map`: optional `Map` applied to `(coeffs_cell, offsets_cell)` before slave selection
- `on_conflict`, `T`: same as the cell-array constructor
"""
function ConstraintHandler(
  space::FESpace,
  space_test::FESpace,
  lhs::Function,
  rhs::Function;
  local_map   = nothing,
  on_conflict = nothing,
  T           = Float64
)
  u     = get_trial_fe_basis(space)
  v     = get_fe_basis(space_test)
  lhs_c = lhs(u, v)
  rhs_c = rhs(v)

  @check CellData.is_matrix_contribution(lhs_c) && CellData.is_vector_contribution(rhs_c) "Expected lhs(u,v) → matrix contribution and rhs(v) → vector contribution"
  @check num_domains(lhs_c) == num_domains(rhs_c) == 1 "Only single-domain contributions are supported"
  trian = first(get_domains(lhs_c))
  @check first(get_domains(rhs_c)) === trian "lhs and rhs must share the same triangulation"

  cell_coeffs  = get_contribution(lhs_c, trian)
  cell_offsets = get_contribution(rhs_c, trian)
  cell_dof_ids = get_cell_dof_ids(space, trian)
  cell_nc      = Arrays.lazy_collect(lazy_map(length, get_cell_dof_ids(space_test, trian)))

  return _build_from_cells(space, cell_dof_ids, cell_coeffs, cell_offsets, cell_nc, T; local_map, on_conflict)
end

# ---------------------------------------------------------------------------
# Constructor 2 — cell, linear: lhs(u), rhs()
# ---------------------------------------------------------------------------

"""
    ConstraintHandler(space, lhs, rhs; local_map, on_conflict, T)

Build a `ConstraintHandler` from a linear form `lhs(u)` and scalar functional `rhs()`.
One slave DOF is selected per cell via `select_constraint_slaves`.

# Arguments
- `space`: trial FE space
- `lhs`: `u → DomainContribution` — per-cell coefficient vectors
- `rhs`: `() → DomainContribution` — per-cell offset scalars
- `local_map`, `on_conflict`, `T`: same as the bilinear constructor
"""
function ConstraintHandler(
  space::FESpace,
  lhs::Function,
  rhs::Function;
  local_map   = nothing,
  on_conflict = nothing,
  T           = Float64
)
  v     = get_fe_basis(space)
  lhs_c = lhs(v)
  rhs_c = rhs()

  @check CellData.is_vector_contribution(lhs_c) && CellData.is_scalar_contribution(rhs_c) "Expected lhs(u) → vector contribution and rhs() → scalar contribution"
  @check num_domains(lhs_c) == num_domains(rhs_c) == 1 "Only single-domain contributions are supported"
  trian = first(get_domains(lhs_c))
  @check first(get_domains(rhs_c)) === trian "lhs and rhs must share the same triangulation"

  cell_coeffs  = get_contribution(lhs_c, trian)
  cell_offsets = get_contribution(rhs_c, trian)
  cell_dof_ids = get_cell_dof_ids(space, trian)
  cell_nc      = fill(1, num_cells(trian))

  return _build_from_cells(space, cell_dof_ids, cell_coeffs, cell_offsets, cell_nc, T; local_map, on_conflict)
end

# ---------------------------------------------------------------------------
# Constructor 3 — patch, bilinear: lhs(u,v), rhs(v)
# ---------------------------------------------------------------------------

"""
    ConstraintHandler(ptopo, space, space_test, lhs, rhs; local_map, on_conflict, T)

Patch-assembly constructor with bilinear form. Assembles `lhs(u,v)` as a matrix per
patch and `rhs(v)` as a vector per patch, then auto-selects slave DOFs.

# Arguments
- `ptopo`: `PatchTopology` defining the patch structure
- `space`: trial FE space (column DOFs / master candidates)
- `space_test`: test FE space (row DOFs / one constraint per test DOF per patch)
- `lhs`: `(u,v) → DomainContribution`
- `rhs`: `v → DomainContribution`
- `local_map`: optional `Map` applied to `(coeffs_patch, offsets_patch)` before slave selection
- `on_conflict`, `T`: same as the cell-array constructor
"""
function ConstraintHandler(
  ptopo::PatchTopology,
  space::FESpace,
  space_test::FESpace,
  lhs::Function,
  rhs::Function;
  local_map   = nothing,
  on_conflict = nothing,
  T           = Float64
)
  assem        = PatchAssembler(ptopo, space, space_test)
  cell_coeffs  = assemble_matrix(lhs, assem, space, space_test)
  cell_offsets = assemble_vector(rhs, assem, space_test)
  cell_nc      = Arrays.lazy_collect(lazy_map(length, assem.strategy.patch_rows))

  return _build_from_cells(space, assem.strategy.patch_cols, cell_coeffs, cell_offsets, cell_nc, T; local_map, on_conflict)
end

# ---------------------------------------------------------------------------
# Constructor 4 — patch, linear: lhs(u), rhs()
# ---------------------------------------------------------------------------

"""
    ConstraintHandler(ptopo, space, lhs, rhs; local_map, on_conflict, T)

Patch-assembly constructor with linear form. Assembles `lhs(u)` as a vector per patch
and `rhs()` as a scalar per patch, then auto-selects one slave DOF per patch.

# Arguments
- `ptopo`: `PatchTopology` defining the patch structure
- `space`: trial FE space
- `lhs`: `u → DomainContribution`
- `rhs`: `() → DomainContribution`
- `local_map`, `on_conflict`, `T`: same as the patch bilinear constructor
"""
function ConstraintHandler(
  ptopo::PatchTopology,
  space::FESpace,
  lhs::Function,
  rhs::Function;
  local_map   = nothing,
  on_conflict = nothing,
  T           = Float64
)
  assem        = PatchAssembler(ptopo, space, space)
  cell_coeffs  = assemble_vector(lhs, assem, space)
  cell_offsets = assemble_scalar(ptopo, rhs())
  cell_nc      = fill(1, Geometry.num_patches(ptopo))

  return _build_from_cells(space, assem.strategy.patch_cols, cell_coeffs, cell_offsets, cell_nc, T; local_map, on_conflict)
end
